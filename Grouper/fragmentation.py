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
    empty_graph = ((), (), (), ())
    atomsRange = range(mol.GetNumAtoms())
    matchesList = [((), tuple(atomsRange), copy.copy(empty_graph))]
    for i, query in enumerate(queries):
        # print(f"\nQuerying {[a.GetSymbol() for a in query.GetAtoms()]}")
        updateMatches = set()
        for current_match in matchesList:
            # should return [match] if query doesn't change anything
            # should return list of matches, corresponding to multiple query options
            new_matches = _generate_new_matches(mol, query, current_match, i)
            # print(f"{new_matches=}")
            for new_match in new_matches:
                updateMatches.add(
                    new_match
                )  # hash function handles adding same match twice
        if updateMatches:  # only update if matches have been found
            matchesList = list(
                updateMatches
            )  # reset matchesList for next round of matching

    # print(f"\n{matchesList=}\n{len(matchesList)=}") # TODO: MatchesList is flaky, order is not preserved

    if incompleteGraphHandler == "remove" or incompleteGraphHandler == "raise error":
        fullyMatchedList = []
        for match in matchesList:
            if not match[1]:
                fullyMatchedList.append(match)
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
            match
        ) in matchesList:  # loop to select smallest possible number of nodes matched
            if len(match[2][1]) < n_nodes:
                n_nodes = len(match[2][1])
                matchesList = [match]  # only use best match
    allGroupGraphList = [
        _group_graph_from_match(
            possible_fragmentation, nodeDefs, incompleteGraphHandler
        )
        for possible_fragmentation in matchesList
    ]
    allGroupGraphList = [groupG for groupG in allGroupGraphList if groupG is not None]
    return allGroupGraphList


def _update_matches_by_state(mol, groupsList, match, groupStr, queryIndex):
    """Take a list of groups, check compatibility with current matching state of molecule, and try to add to state."""
    # matched atoms, unmatched atoms, ((grouped_atoms), (grouped_atom_indices), (group_indices), (grouped_bonds))
    updated_match = ((), (), ((), (), (), ()))
    matched_atoms = match[0]
    unmatched_atoms = match[1]
    groups = match[2][0]
    groupAtomIndices = match[2][1]
    groupIndex = match[2][2]
    bondsTuple = match[2][3]
    for group in groupsList:
        # skip group if atoms in group already matched
        groupFlag = False  # flag for catching if group matched already
        for atom in group:
            if atom in match[0]:
                groupFlag = True  # this group has already been identified
        if groupFlag:
            continue
        matched_atoms += group
        unmatched_atoms = tuple(filter(lambda y: y not in group, unmatched_atoms))
        groups += (groupStr,)
        groupAtomIndices += (group,)
        groupIndex += (queryIndex,)
        groupAttachments = _get_group_bonds(mol, group, groupAtomIndices)
        bondsTuple += groupAttachments
    return (
        matched_atoms,
        unmatched_atoms,
        (groups, groupAtomIndices, groupIndex, bondsTuple),
    )  # return a fresh match


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
            )
            if bonded_groupIndex != -1:
                bonds.append(
                    ((newest_groupIndex, i), (bonded_groupIndex, bonded_hubIndex))
                )  # add a new bond
    return tuple(bonds)


def _generate_new_matches(mol, query, match, queryIndex, returnHandler="all matches"):
    # Substructs are all possible gruops
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

    new_matchesList = []
    for state in iterPossibleStates:
        groupStr = tuple(
            a.GetSymbol() for a in query.GetAtoms()
        )  # used for debugging matches
        # print(f"{state=}")
        new_match = _update_matches_by_state(
            mol, state, match, groupStr, queryIndex
        )  # return `match` if no match possible
        new_matchesList.append(new_match)
    return new_matchesList  # list of hashable matches


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


def _group_graph_from_match(query_match, nodesList, incompleteGraphHandler):
    """Take a matched group object and return a GroupGraph."""
    groupG = GroupGraph()
    groups = query_match[2][2]  # change it so this is group indexes
    # ("C", "C"), ("C", "O") : (0,1), (2,3) ((0,1), (1,0))
    # Group from groupTuple -> name, ports, hubs,
    for group in groups:  # group is an index from the original nodes
        group_node = nodesList[group]
        groupG.add_node(
            group_node.type, group_node.pattern, group_node.hubs, group_node.is_smarts
        )
    # Add edges
    edges = query_match[2][3]
    for edge in edges:
        # need port, not hub for edges
        group0 = edge[0][0]  # this info should be saved in bonds tuple
        group1 = edge[1][0]
        port0 = _get_next_available_port_from_hub(groupG, group0, edge[0][1])
        port1 = _get_next_available_port_from_hub(groupG, group1, edge[1][1])
        if port0 is not None and port1 is not None:
            groupG.add_edge((group0, port0), (group1, port1))
        elif incompleteGraphHandler == "raise error":
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
    # raise ValueError(f"In {groupG=}, no available {hub=} for {group=}, for Group {nodeIndex}")


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
