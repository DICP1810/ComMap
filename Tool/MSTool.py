import sys

import numpy as np


def toolGetSortIndex(mass_list):
    temp_index = {}
    counter = 0
    for i, list_item in enumerate(mass_list):
        for j, _ in enumerate(list_item):
            temp_index[counter] = (i, j)
            counter += 1
    temp = np.array([item for list_item in mass_list for item in list_item])
    result = temp[np.argsort(temp)]
    result_index = [temp_index[order] for order in np.argsort(temp)]
    return result, result_index


def toolFilteredByMass(raw_list, low_bound, up_bound):
    raw_list = np.array(raw_list)
    return np.where((raw_list < low_bound) | (raw_list > up_bound))[0]


def toolInsertDictionary(dictionary, new_key, new_value):
    if new_key in dictionary:
        dictionary[new_key] += new_value
    else:
        dictionary[new_key] = new_value


def toolInsertListItemDictionary(dictionary, new_key, new_value):
    if new_key in dictionary:
        dictionary[new_key].append(new_value)
    else:
        dictionary[new_key] = [new_value]


def toolDecompositionFormular(formular):
    result = [[], []]
    num_str = ''
    for item in formular:
        if item.isalpha():
            result[0].append(item)
        elif item == ')':
            result[1].append(int(num_str))
            num_str = ''
        elif item != '(':
            num_str += item
    return result


def toolGetWord1(inputString, d1, d2):
    start = 0
    end = len(inputString)

    for i in range(len(inputString)):

        if inputString[i] == d1:
            start = i + 1

        if inputString[i] == d2:
            end = i

    return inputString[start:end]


def toolGetWord(input_string, index, d):
    if d is not None and input_string[0] != d:
        input_string = d + input_string
    if d is not None and input_string[-1] != d:
        input_string = input_string + d
    return input_string.split(d)[index + 1]


def toolStr2List(input_string, input_separator):
    if input_string[-1] != input_separator:
        input_string = input_string + input_separator
    return [float(item.strip()) for item in input_string.split(input_separator) if item.strip() != '']


def toolUpdateModificationDictionary(dictionary1, dictionary2):
    for item in dictionary2:
        if item in dictionary1:
            dictionary1[item] += dictionary2[item]
        else:
            dictionary1[item] = dictionary2[item]


def toolAddModificationDictionary(dictionary1, dictionary2):
    new_dictionary = dictionary1.copy()
    toolUpdateModificationDictionary(new_dictionary, dictionary2)
    return new_dictionary


def toolCxDecoyType(protein_list1, protein_list2):
    protein1_decoy = False
    protein2_decoy = False
    for protein1 in protein_list1:
        if 'DECOY_' in protein1:
            protein1_decoy = True
            break
    for protein2 in protein_list2:
        if 'DECOY_' in protein2:
            protein2_decoy = True
            break
    if protein1_decoy and protein2_decoy:
        return 'DD'
    elif protein1_decoy or protein2_decoy:
        return 'TD'
    else:
        return 'TT'
