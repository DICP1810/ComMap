import sys

import numpy as np


def toolGetSortIndex(mass_list):
    """
    对锯齿状二维数据进行排序，获取编号
    :param mass_list: 质量列表
    :return: 排好序的质量列表，索引列表（从大到小）
    """
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
    """
    将原始列表中，返回不处于low和up范围内的值的索引
    :param raw_list: 原始列表
    :param low_bound: 过滤的下界
    :param up_bound: 过滤的上界
    :return: 索引列表，每一项为不符合条件的raw_list索引值
    """
    raw_list = np.array(raw_list)
    return np.where((raw_list < low_bound) | (raw_list > up_bound))[0]


def toolInsertDictionary(dictionary, new_key, new_value):
    """
    将新值插入到字典中，如果新的键值已存在于字典当中，则将新值累加到原来的旧值上
    :param dictionary: 字典，字典的值必须为可以进行‘+’操作的对象
    :param new_key: 要插入的数据的键值
    :param new_value: 新的值
    :return: 更新过的字典
    """
    if new_key in dictionary:
        dictionary[new_key] += new_value
    else:
        dictionary[new_key] = new_value


def toolInsertListItemDictionary(dictionary, new_key, new_value):
    """
    将新值插入到字典中，如果新的键值已存在于字典当中，则将新值累加到原来的旧值的列表中
    :param dictionary: 字典，字典的值必须为列表
    :param new_key: 要插入的数据的键值
    :param new_value: 新的值
    :return: 更新过的字典
    """
    if new_key in dictionary:
        dictionary[new_key].append(new_value)
    else:
        dictionary[new_key] = [new_value]


def toolDecompositionFormular(formular):
    """
    分解分子结构式
    :return:一个嵌套列表，嵌套列表的第一个列表为元素名称，第二个列表为元素的个数
    """
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
    """
    提取出以d作为分隔符的第index项目的值
    :param input_string:待分割的字符串
    :param index:分割后的索引
    :param d:分隔符
    :return:以d作为分割符的第index项
    """
    if d is not None and input_string[0] != d:
        input_string = d + input_string
    if d is not None and input_string[-1] != d:
        input_string = input_string + d
    return input_string.split(d)[index + 1]


def toolStr2List(input_string, input_separator):
    """
    讲输入字符串以分隔符进行分割，返回分割后的每一项，并将其转换为浮点数
    :param input_string: 输入的字符串
    :param input_separator: 输入的分隔符
    :return: 浮点数列表
    """
    if input_string[-1] != input_separator:
        input_string = input_string + input_separator
    return [float(item.strip()) for item in input_string.split(input_separator) if item.strip() != '']


def toolUpdateModificationDictionary(dictionary1, dictionary2):
    """
    将2号字典中的修饰覆盖加到1号字典当中，返回新的修饰字典
    :param dictionary1: 1号修饰字典
    :param dictionary2: 2号修饰字典
    :return: None
    """
    for item in dictionary2:
        if item in dictionary1:
            dictionary1[item] += dictionary2[item]
        else:
            dictionary1[item] = dictionary2[item]


def toolAddModificationDictionary(dictionary1, dictionary2):
    """
    将1号和2号字典中的修饰进行合并，返回新的字典
    :param dictionary1: 1号修饰字典
    :param dictionary2: 2号修饰字典
    :return: 新的修饰字典
    """
    new_dictionary = dictionary1.copy()
    toolUpdateModificationDictionary(new_dictionary, dictionary2)
    return new_dictionary


def toolCxDecoyType(protein_list1, protein_list2):
    """
    判断工字型交联肽段的正库-反库类型
    如果两条肽段都是正库的话，返回TT
    如果两条肽段都是反库的话，返回DD
    如果两条肽段有一条属于正库，有一条属于反库，返回TD
    :param protein_list1: 肽段1所属于的蛋白列表
    :param protein_list2: 肽段2所属于的蛋白列表
    :return: TT或TD或DD
    """
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
