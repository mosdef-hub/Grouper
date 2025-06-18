"""Module for breaking down smiles molecules into Groups and GroupGraphs."""

import copy
from collections import Counter, OrderedDict
from typing import List, Literal, Union

from rdkit import Chem

from Grouper import Group, GroupGraph, exhaustive_fragment


def fragment(
    smiles: str,
    nodeDefs: Union[list, set, tuple],
    returnHandler: Literal["ideal", "quick", "exhuastive"] = "ideal",
    nodeDefsSorter: Literal["size", "priority", "list"] = "list",
    incompleteGraphHandler: Literal["remove", "keep", "raise error"] = "remove",
    matchHubs: bool = False,
) -> list:
    """Fragmente a smiles defined molecule based on a list of nodes.

    Parameters
    ----------
    smiles : string
        A full molecule smiles string
    nodeDefs : list of Grouper.Group> or list of Grouper.GroupExtension
        A set of Groups or GroupExtensions to use for fragmenting molecules.
    returnHandler : string, default "ideal"
        How to handle multiple matches. Options are: "ideal", "quick", "exhaustive"
        - "ideal": return a list of matches that are ideal "i.e. the shortes possible matches using the given order of nodeDefs.
            The list that is returned are sorted on the shortest `GroupGraph.Groups`.
        - "quick": return the match that is found first, discarding the rest along the way. This can lead to the possibility of missing a
            possible match as the acceptable match might get discarded along the way.
        - "exhaustive": return all possible matches, which accounts for possible symmetry matches. Is much slower but more robust than the
            other two returnHandlers.
    nodeDefsSorter : string, default "size"
        How to sort the nodeDefs. Options are: "size", "priority"
        - "size": sort by size of the group, where size is defined first as the number
        of heavy atoms, then the total molecular weight of the group.
        - "priority": sort by priority of the group, where priority is defined by the groupExtension.priority
        attribute. This can be set by the user. If two groups have the same priority, their index in nodeDefs is used.
        - "list": sort by the order of the groups in the list.
    incompleteGraphHandler : string, default "remove"
        How to handle incomplete graphs. Options are: "remove", "keep", "raise error"
        - "remove": remove incomplete graphs, aka graphs that are not fully connected or fragmented.
        - "keep": keep incomplete graphs
        - "raise error": raise an error if no complete graphs are found
    matchHubs : bool, default False
        Whether or not to parse the hubs into the SMARTS string of the Group in nodeDefs is an iterable of Grouper.Groups.
        By default, only the SMARTS or SMILES string of the Group is used for RDKit substructure matching. To make matches also
        take into account where hubs have been positioned around the Group, set this to True.
        Note: This will only work if the Group is an iterable of Grouper.Groups, and the Groups are written as SMARTS.
        Note: This will not be accounted for if returnHandler is set to "exhaustive", which automatically internally
            accounts for the hubs connectivity.


    Notes
    -----
    1. `matchHubs` can result in different matching behavior, since there are not enough
        hubs to for the middle carbon to form two bonds in the example below.
    ```python
    group1 = Grouper.Group("group1", "[C]", [0], True)
    fragment("CCC", [group1], matchHubs=True) == []
    gG = Grouper.GroupGraph()
    gG.add_node("group1", "[C]", [0]*4, True)
    gG.add_node("group1", "[C]", [0]*4, True)
    gG.add_node("group1", "[C]", [0]*4, True)
    gG.add_edge((0,0), (1,0), 1)
    gG.add_edge((2,0), (1,1), 1)
    fragment("CCC", [group1], matchHubs=False) == GroupGraph
    ```

    2. returnHandler="exhaustive" can be slow for large matching structures due to accounting for all
    possible symmetries of the groups supplied in `nodeDefs`.

    3. If a list of SMARTS is supplied to `nodeDefs` instead of a list of Groups, then the final graphs
    will have nodes where every possible hub is added to the Group based on the subatoms valency and connectivity.
    This is recommended for the most intuitive matching behavior.

    Returns
    -------
    list of Grouper.GroupGraph
        A list of GroupGraphs that represent the possible fragmentations of the molecule.
    """
    # some preliminary type checking and error handling
    # Handle both GroupExtension and Groups
    mol = Chem.MolFromSmiles(smiles)  # convert to RDKit mol
    if not smiles or not nodeDefs:  #
        return []
    if returnHandler == "exhaustive":  # call cpp fragmenter for exhaustive matching
        _, nodeDefs = _generate_queries_from_nodedefs(
            nodeDefs, nodeDefsSorter, matchHubs
        )
        return exhaustive_fragment(smiles, set(nodeDefs))
    queries, nodeDefs = _generate_queries_from_nodedefs(
        nodeDefs, nodeDefsSorter, matchHubs
    )

    matchesList = [MatchState(mol.GetNumAtoms())]
    for i, query in enumerate(queries):
        updateMatchesDict = OrderedDict()  # use dict to keep order of set
        for currentMatchState in matchesList:
            # should return [match] if query doesn't change anything
            # should return list of matches, corresponding to multiple query options
            newMatchesList = _generate_new_matches(
                mol, query, currentMatchState, i, returnHandler
            )
            for newMatchState in newMatchesList:
                updateMatchesDict.update(
                    {newMatchState: None}
                )  # matchState is hashable
        if updateMatchesDict:  # only update if matches have been found
            matchesList = list(
                updateMatchesDict.keys()
            )  # reset matchesList for next round of matching
    # Filter all matches based on fragment arguments
    if incompleteGraphHandler == "remove" or incompleteGraphHandler == "raise error":
        fullyMatchedList = []
        for matchState in matchesList:
            if matchState.isFullyFragmented():
                fullyMatchedList.append(matchState)
        matchesList = fullyMatchedList
        if incompleteGraphHandler == "raise error" and not matchesList:
            raise ValueError(
                "No complete graphs found. Please use a different set of nodDefs."
            )
    # sort by least number of nodes
    matchesList.sort(key=lambda x: len(x.all_group_indices))
    allGroupGraphList = [
        _group_graph_from_match(
            possible_fragmentation, nodeDefs, incompleteGraphHandler
        )
        for possible_fragmentation in matchesList
    ]
    allGroupGraphList = [groupG for groupG in allGroupGraphList if groupG is not None]
    return allGroupGraphList


class MatchState:
    __slots__ = (
        "all_matched_atoms",
        "all_unmatched_atoms",
        "all_group_names",
        "all_group_indices",
        "nodeDefs_indices",
        "group_bonds",
    )

    def __init__(self, n_atoms):
        self.all_matched_atoms = ()
        self.all_unmatched_atoms = tuple(range(n_atoms))
        self.all_group_names = ()
        self.all_group_indices = ()
        self.nodeDefs_indices = ()
        self.group_bonds = ()

    def isFullyFragmented(self):
        return not self.all_unmatched_atoms

    def __eq__(self, other):
        return (
            self.all_matched_atoms == other.all_matched_atoms
            and self.all_unmatched_atoms == other.all_unmatched_atoms
            and self.all_group_names == other.all_group_names
            and self.all_group_indices == other.all_group_indices
            and self.nodeDefs_indices == other.nodeDefs_indices
            and self.group_bonds == other.group_bonds
        )

    def __hash__(self):
        return hash(
            (
                self.all_matched_atoms,
                self.all_unmatched_atoms,
                self.all_group_names,
                self.all_group_indices,
                self.nodeDefs_indices,
                self.group_bonds,
            )
        )

    def __repr__(self):
        return (
            f"{self.all_unmatched_atoms=}\n{self.all_group_names=}\n{self.group_bonds=}"
        )


def _update_matches_by_state(mol, groupsList, oldMatchState, groupStr, nodeDefsIndex):
    matchState = copy.deepcopy(oldMatchState)
    """Take a list of groups, check compatibility with current matchState molecule, and update slots."""
    for group in groupsList:
        # skip group only if atoms in group already used in another group
        groupFlag = False  # flag for catching if group matched already
        for atom in group:
            if atom in matchState.all_matched_atoms:
                groupFlag = True  # this group has already been used
        if groupFlag:
            continue  # at least one previous atom was already found
        matchState.all_matched_atoms += group
        matchState.all_unmatched_atoms = tuple(
            filter(lambda y: y not in group, matchState.all_unmatched_atoms)
        )
        matchState.all_group_names += (groupStr,)
        matchState.all_group_indices += (group,)
        matchState.nodeDefs_indices += (nodeDefsIndex,)
        matchState.group_bonds += _get_group_bonds(
            mol, group, matchState.all_group_indices
        )
    return matchState  # return a fresh match


def _identify_bond_partner_group_and_hub(currentGroups, bondedAtomIndex):
    """Return -1, -1 if no bond to be made"""

    for j, group_by_index in enumerate(currentGroups):
        for k, atomid in enumerate(group_by_index):
            if atomid == bondedAtomIndex:
                return (j, k)
    return -1, -1


def _get_group_bonds(mol, group, latest_group_indices):
    bonds = []  # return list
    newest_groupIndex = (
        len(latest_group_indices) - 1
    )  # updating based on index of last group. Issue could be if two of the same groups are added at once
    for i, atomIndex in enumerate(group):
        atom = mol.GetAtomWithIdx(atomIndex)
        # ((groupIndex1, port1Hub), (groupIndex2, group2Hub))
        for bond in atom.GetBonds():
            if atomIndex == bond.GetBeginAtomIdx():  # query is begin atom
                bondedAtomIndex = bond.GetEndAtomIdx()
            else:  # query is end atom
                bondedAtomIndex = bond.GetBeginAtomIdx()
            # take bondedAtom, find group Index and group Hub
            bonded_groupIndex, bonded_hubIndex = _identify_bond_partner_group_and_hub(
                latest_group_indices[:-1], bondedAtomIndex
            )  # latest_group_indices[:-1] because we don't want to count group, which is already added to latest_group_indices
            if bonded_groupIndex != -1:
                bond_order = int(
                    bond.GetBondTypeAsDouble()
                )  # NOTE: Since we have ints, aromatic 1.5 bonds are treated as single bonds
                bonds.append(
                    (
                        (newest_groupIndex, i),
                        (bonded_groupIndex, bonded_hubIndex),
                        bond_order,
                    )
                )  # add a new bond
    return tuple(bonds)


def _generate_new_matches(mol, query, matchState, nodeDefsIndex, returnHandler="ideal"):
    # Substructs are all possible groups
    substructsList = mol.GetSubstructMatches(query)
    # Generate all maximal sets of compatible groups
    # i.e. [(0,1), (0,3)] -> [[(0,1)], [(0,3)]], since both have index 0.
    # list of states, which is a list of substructures
    if returnHandler.lower() == "quick":
        iterPossibleStates = _get_first_compatible_tuples(substructsList)
    elif returnHandler.lower() == "ideal":
        iterPossibleStates = _get_maximal_compatible_tuples(substructsList)
    else:
        raise ValueError(
            f"Invalid returnHandler: {returnHandler}. Must be one of `ideal`, `quick`, or `exhaustive`"
        )

    newMatchesList = []
    for state in iterPossibleStates:
        groupStr = tuple(
            a.GetSymbol() for a in query.GetAtoms()
        )  # The atom strings that will be used to define the group
        newMatch = _update_matches_by_state(
            mol, state, matchState, groupStr, nodeDefsIndex
        )  # return `match` if no match possible
        newMatchesList.append(newMatch)
    return newMatchesList  # list of hashable matches


def _get_maximal_compatible_tuples(tupleofGroups):
    """
    Only return longest possible compatible structures. Return [] if no compatible Groups
    """
    if not tupleofGroups:
        return []
    compatibleGroupsList = [(tupleofGroups[0],)]
    for group in tupleofGroups[1:]:  # iter over remaining groups
        loopUpdateGroupsSet = set()
        for groupTuple in compatibleGroupsList:  # iter over previously paired groups
            # check compatibility function
            loopUpdateGroupsSet.add(groupTuple)
            loopUpdateGroupsSet.add((group,))
            clashes = _mask_number_clashes(
                groupTuple, group
            )  # indexList that corresponds to clashes
            compatibleGroups = tuple(
                [group for group, mask in zip(groupTuple, clashes) if not mask]
            )
            loopUpdateGroupsSet.add(compatibleGroups + (group,))

        # take longest element of set
        compatibleGroupsList = list(loopUpdateGroupsSet)
    maxLength = max(len(compatibleGroups) for compatibleGroups in compatibleGroupsList)
    compatibleGroupsList = list(
        compatibleGroups
        for compatibleGroups in compatibleGroupsList
        if len(compatibleGroups) == maxLength
    )  # keep only longest options
    return compatibleGroupsList


def _get_first_compatible_tuples(tupleofGroups):
    """
    Only return first matched structures. Return [] if no compatible Groups
    """
    if not tupleofGroups:
        return []
    # compatibleGroupsList = [(tupleofGroups[0],)]
    compatibleGroupsList = [tupleofGroups[0]]
    for group in tupleofGroups[1:]:  # iter over remaining groups
        # check compatibility function
        clashes = _mask_number_clashes(
            compatibleGroupsList, group
        )  # indexList that corresponds to clashes
        n_clashes = sum(clashes)
        if not n_clashes:
            compatibleGroupsList.append(group)

    return [tuple(compatibleGroupsList)]  # keep in general format


def _get_exhaustive_compatible_tuples(tupleofGroups):
    """
    This function is currently only partially exhaustive.
    Please rewrite to be more inline with getMaximalCompatibleTuples.
    For fully testing, try getExhaustiveCompatibleSubSets(((0,1),(1,2),(2,3)))
    which returns (0,1),(2,3) : (1,2) : (2,3), but not (0,1) or ()

    """
    # sort
    # initialize output datastructure -> [[group1, group2], [group2, group3]]

    compatibleGroupsList = [[tupleofGroups[0]]]
    for group in tupleofGroups[1:]:  # iter over remaining groups
        loopUpdateGroupsList = []
        for groupList in compatibleGroupsList:  # iter over previously paired groups
            # check compatibility function
            clashes = _mask_number_clashes(
                groupList, group
            )  # indexList that corresponds to clashes
            n_clashes = sum(clashes)

            if n_clashes > 1:  # too many clashes, don't combine
                loopUpdateGroupsList.append(groupList)  # no changes
                if [group] not in loopUpdateGroupsList:
                    loopUpdateGroupsList.append([group])
            elif n_clashes == 1:
                loopUpdateGroupsList.append(groupList)  # keep old
                newList = [elem for elem, mask in zip(groupList, clashes) if not mask]
                if not newList:  # handle repeat adding just a single group
                    # This is touchy. Only add a single group if all other possible matches are 1 in length.
                    if [group] not in loopUpdateGroupsList:
                        loopUpdateGroupsList.append([group])
                else:
                    loopUpdateGroupsList.append(
                        newList + [group]
                    )  # keep mask plus new group
            elif not n_clashes:  # add group
                loopUpdateGroupsList.append(groupList + [group])
        # need to update compatibleGroupsList for next cycle
        compatibleGroupsList = loopUpdateGroupsList

    compatibleGroupsList = max(
        compatibleGroupsList, key=len
    )  # keep all longest options

    return compatibleGroupsList


def _mask_number_clashes(groupList, group):
    """Return a mask of 0's and 1's if there are clashes with group."""
    groupSet = set(group)
    maskList = []
    for compGroup in groupList:
        compSet = set(compGroup)
        maskList.append(not compSet.isdisjoint(groupSet))
    return maskList


def _group_graph_from_match(matchState, nodesList, incompleteGraphHandler):
    """Take a MatchState object and return a GroupGraph."""
    groupG = GroupGraph()
    # Group from matchState.nodeDefs_indices, which tells you which groups were found
    groupsTuple = matchState.nodeDefs_indices
    for group in groupsTuple:  # group is an index from the original nodes
        group_node = nodesList[group]
        groupG.add_node(
            group_node.type,
            group_node.pattern,
            group_node.hubs,
            group_node.pattern_type,
        )
    # Add edges
    edgesTuple = matchState.group_bonds
    for (group0, hub0), (group1, hub1), bond_order in edgesTuple:
        # need the port, but have the hub for creation of GroupGraph edges
        port0 = _get_next_available_port_from_hub(groupG, group0, hub0)
        port1 = _get_next_available_port_from_hub(groupG, group1, hub1)
        if port0 is not None and port1 is not None:
            groupG.add_edge((group0, port0), (group1, port1), bond_order)
        elif incompleteGraphHandler == "raise error":
            # raise errors if we generated a failing graph
            if port0 is None:
                raise ValueError(
                    f"No more ports in {group0=} to connect to {group1=} for {groupG=}"
                )
            else:
                raise ValueError(
                    f"No more ports in {group1=} to connect to {group0=} for {groupG=}"
                )
        elif incompleteGraphHandler == "remove":
            return None  # just skip graphs with missing edges
    return groupG


def _get_next_available_port_from_hub(groupG, nodeIndex, hub):
    used_ports = []
    for edge in groupG.edges:  # do we handle double bonds
        if edge[0] == nodeIndex:
            used_ports.append(edge[1])
        if edge[2] == nodeIndex:
            used_ports.append(edge[3])

    group = groupG.nodes[nodeIndex]
    for portIndex, checkHub in enumerate(group.hubs):
        if hub == checkHub and portIndex not in used_ports:
            return portIndex


def _get_hubs_from_string(pattern):
    letterIndex = 0
    hubs = []
    charValencyMap = {
        "C": 4,
        "N": 3,
        "O": 2,
        "F": 1,
        "P": 3,
        "B": 3,
        "Cl": 1,
        "I": 1,
        "Li": 1,
        "Si": 4,
        "S": 2,
        "Se": 2,
        "Br": 1,
        "H": 1,
    }

    # main working loop
    # subtract 1 from charValencyMap since only looking for hubs to potential groups outside of the group
    # at minimum there is 1 bond within the group
    for char1, char2 in zip(pattern, pattern[1:]):
        if (
            char1 + char2 in charValencyMap
        ):  # matched a pair of two letters, so don't match on just a single letter
            hubs.extend([letterIndex for _ in range(charValencyMap[char1 + char2] - 1)])
            letterIndex += 1
        elif char1 in charValencyMap:  # matched a single letter
            hubs.extend([letterIndex for _ in range(charValencyMap[char1] - 1)])
            letterIndex += 1
    # handle last string
    if pattern[-1] in charValencyMap:
        if pattern[-1] in charValencyMap:
            hubs.extend([letterIndex for _ in range(charValencyMap[pattern[-1]] - 1)])
            letterIndex += 1

    # handle edge case where only 1 atom in pattern, so no -1 in pattern:
    if letterIndex == 1:
        hubs.append(0)  # add in one more possible bonding state

    return hubs


def _smarts_with_ports(smarts: str, hubs: List[int]):
    """Take a SMARTS string, and return the SMARTS that match the specified hubs in the group.
    i.e smarts="[C]", hubs=[0] is a CH3 or CH4, so the explicit_smarts should be "[CH4,CH3]"
    as opposed to ["C"], hubs=[0,0,0], which is explicit_smarts "[CH1,CH2,CH3,CH4]"
    Note: "[CH1]" will not match CH2 and CH3 queries
    """
    standardElementValencyMap = {
        "H": 1,
        "B": 3,
        "C": 4,
        "N": 3,
        "O": 2,
        "F": 1,
        "P": 3,
        "S": 2,
        "Cl": 1,
        "Br": 1,
        "I": 1,
    }
    mol = Chem.MolFromSmarts(smarts)
    num_portsCounter = Counter(hubs)
    new_smartsStr = ""
    smarts_atomsList = _split_smarts_by_atom(smarts)  # each atom's smarts string
    for i, (atom, atomSmarts) in enumerate(zip(mol.GetAtoms(), smarts_atomsList)):
        atomStr = atomSmarts
        replaceStr = ""
        symbol = atom.GetSymbol()
        atomStr = _handle_branching(atomStr)  # handle branching
        atomStr = _handle_and_statements(
            atomStr, symbol
        )  # put and statements first in []
        valence = standardElementValencyMap.get(atom.GetSymbol(), None)
        bonds_occupied = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        n_ports = num_portsCounter[i]
        # 2 hubs-> CH2, CH1, CH0 # ports means at least valence - bonds - ports hydrogens
        # 1 hubs -> CH2, CH1
        # 0 hubs -> CH2
        maxn_H = int(
            valence - bonds_occupied
        )  # max hydrogens is valence - bonds_occupied, min is valence - bonds - n_ports
        for j in range(maxn_H - n_ports, maxn_H + 1):
            replaceStr += symbol + "H" + str(j) + ","
        replaceStr = replaceStr[:-1]  # remove last comma
        atomStr = atomStr.replace(symbol, replaceStr)
        new_smartsStr += atomStr
    return new_smartsStr


def _handle_branching(atomStr):
    """Handle branching in a smarts string."""
    innermost_left_paren = atomStr.rfind("(")
    if innermost_left_paren == -1:
        outStr = "["
    else:
        outStr = atomStr[: innermost_left_paren + 1] + "["
    innermost_right_paren = atomStr.find(")")
    if innermost_left_paren == -1:
        outStr += atomStr[innermost_left_paren + 1 :] + "]"
    else:
        outStr += (
            atomStr[innermost_left_paren + 1 : innermost_right_paren]
            + "]"
            + atomStr[innermost_right_paren:]
        )
    return outStr


def _handle_and_statements(smarts, symbol):
    """Will only handle ;, not & or ,"""
    indices_left = smarts.count("(")
    indices_right = smarts.count(")")

    splitSmarts = smarts.strip("()[]").split(";")
    translation_table = dict.fromkeys(map(ord, "[]()"), None)
    outStr = "["
    for subString in splitSmarts:
        if symbol in subString:
            continue
        else:
            outStr += subString.translate(translation_table) + ";"
    outStr += symbol + "]"  # add element symbol last
    # add back in parantheses
    outStr = indices_left * "(" + outStr + indices_right * ")"
    return outStr


def _split_smarts_by_atom(smarts):
    """Split a smarts string into a list for each atom."""
    standardElementValencyMap = {
        "H": 1,
        "B": 3,
        "C": 4,
        "N": 3,
        "O": 2,
        "F": 1,
        "P": 3,
        "S": 2,
        "Cl": 1,
        "Br": 1,
        "I": 1,
    }
    atomsList = [""]
    added_atom_once = False
    cleanSmarts = smarts.replace("[", "").replace("]", "")
    for char in cleanSmarts:
        if char == "(":
            added_atom_once = False
            atomsList.append("(")
            continue
        elif char == ")":
            atomsList[-1] += ")"
            atomsList.append("")
            added_atom_once = False
            continue
        elif char in standardElementValencyMap and not added_atom_once:
            # only add once, stay here
            added_atom_once = True
        elif char in standardElementValencyMap and added_atom_once:
            atomsList.append(char)
            continue
        atomsList[-1] += char

    # handle extra bits
    outSymbolsList = []
    i = 0
    while i < len(atomsList):
        symbol = atomsList[i]
        if symbol == "":
            pass
        elif symbol == ")":
            outSymbolsList[-1] += symbol
        elif symbol == "(":
            outSymbolsList.append(symbol)
            i += 1
            outSymbolsList[-1] += atomsList[i]
        else:
            outSymbolsList.append(symbol)
        i += 1

    return outSymbolsList


def _generate_queries_from_nodedefs(nodeDefs, nodeDefsSorter, matchHubs):
    """Generate a list of RDKit.Chem.Mol queries from the nodeDefs."""
    queries = []
    if not isinstance(nodeDefs, list):  # indexible list
        nodeDefs = list(nodeDefs)
    if isinstance(nodeDefs[0], str):  # create groups
        newNodeDefs = []
        for i, pattern in enumerate(nodeDefs):
            query = Chem.MolFromSmarts(pattern)
            pattern_type = "SMARTS"
            if not query:  # failed at parsing smarts, try smiles
                query = Chem.MolFromSmiles(pattern)
                pattern_type = "SMILES"
            hubs = _get_hubs_from_string(pattern)
            try:
                newNodeDefs.append(Group(f"type {i}", pattern, hubs, pattern_type))
            except ValueError:  # TODO: add better handling of loading patterns to group
                newNodeDefs.append(Group(f"type {i}", pattern, hubs, pattern_type))
            queries.append(query)
        nodeDefs = newNodeDefs
    elif matchHubs:
        for group in nodeDefs:
            if group.pattern_type.upper() == "SMARTS":
                queries.append(
                    Chem.MolFromSmarts(_smarts_with_ports(group.pattern, group.hubs))
                )
            elif group.pattern_type.upper() == "SMILES":
                queries.append(  # only apply hubs if a SMARTS string was given
                    Chem.MolFromSmiles(group.pattern)
                )
            else:
                raise ValueError(
                    f"Group {group} is not a SMILES or SMARTS pattern. Please set `group.pattern_type` to one of these options to use for fragmentation."
                )
    else:  # by defaults just use the smarts given
        for group in nodeDefs:
            if group.pattern_type.upper() == "SMARTS":
                queries.append(Chem.MolFromSmarts(group.pattern))
            elif group.pattern_type.upper() == "SMILES":
                queries.append(Chem.MolFromSmiles(group.pattern))
            else:
                raise ValueError(
                    f"Group {group} is not a SMILES or SMARTS pattern. Please set `group.pattern_type` to one of these options to use for fragmentation."
                )
    # sort queries
    if nodeDefsSorter == "size":
        queries, nodeDefs = (
            list(t)
            for t in zip(
                *sorted(zip(queries, nodeDefs), key=_size_sorter, reverse=True)
            )
        )
    return queries, nodeDefs


def _size_sorter(pair):
    return (
        pair[0].GetNumHeavyAtoms(),
        sum(atom.GetMass() for atom in pair[0].GetAtoms()),
    )
