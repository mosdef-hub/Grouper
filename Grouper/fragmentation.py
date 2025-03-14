"""Module for breaking down smiles molecules into Groups and GroupGraphs."""

import copy
from typing import Literal, Union

from rdkit import Chem

from Grouper import Group, GroupGraph


def fragment(
    smiles: str,
    nodeDefs: Union[list, set, tuple],
    returnHandler: Literal[
        "single match", "all matches", "fast match"
    ] = "single match",
    nodeDefsSorter: Literal["size", "priority", "list"] = "list",
    incompleteGraphHandler: Literal["remove", "keep", "raise error"] = "remove",
) -> list:
    """Fragmenter a smiles definied molecule based on a list of nodes.

    Parameters
    ----------
    smiles : string
        A full molecule smiles string
    nodeDefs : list of Grouper.Group> or list of Grouper.GroupExtension
        A set of Groups or GroupExtensions to use for fragmenting molecules.
    returnHandler : string, default "single match"
        How to handle multiple matches. Options are: "single match", "all matches", "fast match"
        - "single match": return a single match, which is the best match found using the least number of nodes.
        - "fast match": return a single match, which is the best match found without trying all branching of matches.
        - "all matches": return all matches, which is a list of all possible matches.
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

    -------
    list of Grouper.GroupGraph
        A list of GroupGraphs that represent the possible fragmentations of the molecule.
    """
    # some preliminary type checking and error handling
    # Handle both GroupExtension and Groups
    mol = Chem.MolFromSmiles(smiles)  # convert to RDKit mol
    queries = []
    if not smiles or not nodeDefs:  #
        return []
    if not isinstance(nodeDefs, list):  # indexible list
        nodeDefs = list(nodeDefs)
    if isinstance(nodeDefs, list) and isinstance(nodeDefs[0], str):  # create groups
        newNodeDefs = []
        for i, pattern in enumerate(nodeDefs):
            query = Chem.MolFromSmarts(pattern)
            is_smarts = True
            if not query:  # failed at parsing smarts, try smiles
                query = Chem.MolFromSmiles(pattern)
                is_smarts = False
            hubs = _get_hubs_from_string(pattern)
            try:
                newNodeDefs.append(Group(f"type {i}", pattern, hubs, is_smarts))
            except ValueError:  # TODO: add better handling of loading patterns to group
                newNodeDefs.append(Group(f"type {i}", pattern, hubs, is_smarts))
            queries.append(query)
        nodeDefs = newNodeDefs
    else:
        for group in nodeDefs:
            if group.is_smarts:
                queries.append(Chem.MolFromSmarts(group.pattern))
            else:
                queries.append(Chem.MolFromSmiles(group.pattern))
    # sort queries
    if nodeDefsSorter == "size":
        size_sorter = lambda pair: (
            pair[0].GetNumHeavyAtoms(),
            sum(atom.GetMass() for atom in pair[0].GetAtoms()),
        )
        queries, nodeDefs = (
            list(t)
            for t in zip(*sorted(zip(queries, nodeDefs), key=size_sorter, reverse=True))
        )

    matchesList = [MatchState(mol.GetNumAtoms())]
    for i, query in enumerate(queries):
        # print(f"\nQuerying {[a.GetSymbol() for a in query.GetAtoms()]}")
        updateMatchesSet = set()
        print(f"{matchesList=}")
        for currentMatchState in matchesList:
            # should return [match] if query doesn't change anything
            # should return list of matches, corresponding to multiple query options
            newMatchesList = _generate_new_matches(mol, query, currentMatchState, i)
            print(f"{newMatchesList=}")
            for newMatchState in newMatchesList:
                updateMatchesSet.add(newMatchState)  # matchState is hashable
        if updateMatchesSet:  # only update if matches have been found
            matchesList = list(
                updateMatchesSet
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
    if returnHandler == "single match":
        if not matchesList:  # return null if empty
            return []
        n_nodes = mol.GetNumAtoms()
        for (
            matchState
        ) in matchesList:  # loop to select smallest possible number of nodes matched
            if len(matchState.all_group_indices) < n_nodes:
                n_nodes = len(matchState.all_group_indices)
                matchesList = [matchState]  # only use best match
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
                bonded_groupIndex = j
                bonded_hubIndex = k
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


def _generate_new_matches(
    mol, query, matchState, nodeDefsIndex, returnHandler="all matches"
):
    # Substructs are all possible groups
    substructsList = mol.GetSubstructMatches(query)
    # Generate all maximal sets of compatible groups
    # i.e. [(0,1), (0,3)] -> [[(0,1)], [(0,3)]]
    # list of states, which is a list of substructures
    # print(f"{substructsList=}")
    if returnHandler == "fast match":
        iterPossibleStates = _get_first_compatible_tuples(substructsList)
    elif (
        returnHandler == "all matches"
    ):  # may also do getExhaustiveCompatibleSubSets for returnHandler
        iterPossibleStates = _get_maximal_compatible_tuples(substructsList)
    elif returnHandler == "single match":
        iterPossibleStates = _get_maximal_compatible_tuples(substructsList)
    else:
        raise ValueError(
            f"Invalid returnHandler: {returnHandler}. Must be one of `all matches`, `fast match`, or `single match`"
        )

    print(f"##########{iterPossibleStates=}\n\n")
    newMatchesList = []
    for state in iterPossibleStates:
        groupStr = tuple(
            a.GetSymbol() for a in query.GetAtoms()
        )  # used for debugging matches
        # print(f"{state=}")
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
        # print(f"{compatibleGroupsList=}")
        for groupTuple in compatibleGroupsList:  # iter over previously paired groups
            # check compatibility function
            loopUpdateGroupsSet.add(groupTuple)
            loopUpdateGroupsSet.add((group,))
            clashes = _mask_number_clashes(
                groupTuple, group
            )  # indexList that corresponds to clashes
            n_clashes = sum(clashes)
            # print(f"{clashes=}")
            compatibleGroups = tuple(
                [group for group, mask in zip(groupTuple, clashes) if not mask]
            )
            # print(f"{group=}: {compatibleGroups=}")
            loopUpdateGroupsSet.add(compatibleGroups + (group,))

        # take longest element of set
        # print(f"{loopUpdateGroupsSet=}")
        compatibleGroupsList = list(loopUpdateGroupsSet)
    maxLength = max(len(compatibleGroups) for compatibleGroups in compatibleGroupsList)
    compatibleGroupsList = list(
        compatibleGroups
        for compatibleGroups in compatibleGroupsList
        if len(compatibleGroups) == maxLength
    )  # keep only longest options
    # print(f"{compatibleGroupsList=}")
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
        # print(f"{compatibleGroupsList=}")
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
        # print(f"{compSet=} == {groupSet=}")
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
            group_node.type, group_node.pattern, group_node.hubs, group_node.is_smarts
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
