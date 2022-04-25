import os
import numpy as np
import datetime
from math import cos
from System.MSData import CDataPack
from Function.rankAndOutput import rankAndOutput

class CFunctionQCRank:

    def rank(self,data_package=CDataPack()):
        structure_dictionary,pdb_list,ppi_list=self.__read_result(data_package)
        dis_matrix=self.__generate_matrix(structure_dictionary,pdb_list,ppi_list)
        pdb_ppi_result,features=self.__get_feature(pdb_list,ppi_list,dis_matrix,structure_dictionary,data_package)
        score_dictionary=self.__rescore(pdb_ppi_result,features,data_package)
        self.__rankAndOutPut(pdb_ppi_result,score_dictionary,data_package)

    @staticmethod
    def __read_result(data_package):
        output_path = os.path.join(data_package.my_config.E1_PATH_EXPORT,
                                   f'structure_alignment_result_PDB.{data_package.my_config.C66_OUTPUT_DATE}.csv')
        structure_dictionary = {}
        pdb_list = []
        ppi_list = []
        pdb_name = ''
        with open(output_path, 'r') as f:
            for line in f:
                if 'PDB code' in line:
                    continue
                elif line.strip():
                    if line[0] != ',':
                        pdb_name = line.split(',')[0]
                        pdb_list.append(pdb_name)
                        structure_dictionary[pdb_name] = {}
                    else:
                        ppi = line.split(',')[1]
                        distance = float(line.split(',')[5])
                        if ppi not in structure_dictionary[pdb_name]:
                            structure_dictionary[pdb_name][ppi] = [distance]
                        else:
                            structure_dictionary[pdb_name][ppi].append(distance)
                        ppi_list.append(ppi)
        pdb_list = list(set(pdb_list))
        ppi_list = list(set(ppi_list))
        return structure_dictionary,pdb_list,ppi_list

    @staticmethod
    def __generate_matrix(structure_dictionary,pdb_list,ppi_list):
        distance_matrix = np.zeros((len(pdb_list), len(ppi_list)))
        for i, pdb in enumerate(pdb_list):
            for j, ppi in enumerate(ppi_list):
                if ppi in structure_dictionary[pdb]:
                    distance_matrix[i][j] = min(structure_dictionary[pdb][ppi])
                else:
                    distance_matrix[i][j] = np.nan
        return distance_matrix

    @staticmethod
    def __get_feature(pdb_list,ppi_list,distance_matrix,structure_dictionary,data_package=CDataPack()):
        features={}
        # 特征提取 1 同一PDB在多个PPI下的距离进行条件概率计算
        pdb_distance_distribution = []
        for i, pdb in enumerate(pdb_list):
            temp_array = np.zeros(101)
            for dis in distance_matrix[i]:
                if np.isnan(dis):
                    continue
                if dis > 100:
                    temp_array[-1] += 1
                else:
                    temp_array[int(dis)] += 1
            pdb_distance_distribution.append(temp_array / np.sum(temp_array))
        # 特征提取 2 同一PPI在多个PDB下的距离进行条件概率计算
        ppi_distribution = []
        for i, ppi in enumerate(ppi_list):
            temp_array = np.zeros(101)
            for dis in distance_matrix[:, i]:
                if np.isnan(dis):
                    continue
                if dis > 100:
                    temp_array[-1] += 1
                else:
                    temp_array[int(dis)] += 1
            ppi_distribution.append(temp_array / np.sum(temp_array))
        # 读取结果
        pdb_ppi_result={}
        for index in range(len(data_package.my_CxResult.CFLOW6_PPI)):
            for pdb_name in data_package.my_CxResult.CFLOW6_DISTANCE[index]:
                for item in data_package.my_CxResult.CFLOW6_DISTANCE[index][pdb_name]:
                    if pdb_name in pdb_ppi_result:
                        pdb_ppi_result[pdb_name].append((data_package.my_CxResult.CFLOW6_PPI[index],
                                                                pdb_name, item[1], item[2], item[3], item[4],
                                                                str('dist {},/{}/{}/{}/{}/CA,/{}/{}/{}/{}/CA'.format(
                                                                    data_package.my_CxResult.CFLOW6_PPI[index],
                                                                    pdb_name, item[2], item[7],item[5], pdb_name, item[3], item[8],
                                                                    item[6]))))
                    else:
                        pdb_ppi_result[pdb_name]=[(data_package.my_CxResult.CFLOW6_PPI[index],
                                                                pdb_name, item[1], item[2], item[3], item[4],
                                                                str('dist {},/{}/{}/{}/{}/CA,/{}/{}/{}/{}/CA'.format(
                                                                    data_package.my_CxResult.CFLOW6_PPI[index],
                                                                    pdb_name, item[2],item[7], item[5], pdb_name, item[3], item[8],
                                                                    item[6])))]
        for pdb_item in pdb_ppi_result:
            features[pdb_item]=[]
            for info_item in pdb_ppi_result[pdb_item]:
                ppi_index=data_package.my_CxResult.CFLOW6_PPI.index(info_item[0])
                item_distance=info_item[5]

                if '{:.3f}'.format(item_distance)!=str(min(structure_dictionary[pdb_item][str(info_item[0])])):
                    # 这不是此PDB内的最小值
                    features[pdb_item].append([])
                    continue
                else:
                    min_distance_frequency = 0
                    for item in structure_dictionary[pdb_item][str(info_item[0])]:
                        if int(item) == int(item_distance):
                            min_distance_frequency += 1
                    distance_delta = 1 / (1 + abs(item_distance - data_package.my_config.C65_LINKER_ARM_LENGTH))
                if item_distance > data_package.my_config.C65_LINKER_MAX_ARM_LENGTH:
                    pdb_ppi_probability = pdb_distance_distribution[pdb_list.index(pdb_item)][-1]
                    ppi_probability = ppi_distribution[ppi_list.index(str(info_item[0]))][-1]
                else:
                    pdb_ppi_probability = pdb_distance_distribution[pdb_list.index(pdb_item)][int(item_distance)]
                    ppi_probability = ppi_distribution[ppi_list.index(str(info_item[0]))][int(item_distance)]
                features[pdb_item].append([pdb_ppi_probability, ppi_probability, min_distance_frequency, distance_delta])
                '''
                min_distance_frequency = 0
                for item in structure_dictionary[pdb_item][str(info_item[0])]:
                    if int(item)//10 == int(item_distance)//10:
                        min_distance_frequency += 1
                distance_delta = 1 / (1 + abs(item_distance - data_package.my_config.C65_LINKER_MAX_ARM_LENGTH))
                if item_distance > data_package.my_config.C65_LINKER_MAX_ARM_LENGTH:
                    pdb_ppi_probability = pdb_distance_distribution[pdb_list.index(pdb_item)][-1]
                    ppi_probability = ppi_distribution[ppi_list.index(str(info_item[0]))][-1]
                else:
                    pdb_ppi_probability = pdb_distance_distribution[pdb_list.index(pdb_item)][int(item_distance)]
                    ppi_probability = ppi_distribution[ppi_list.index(str(info_item[0]))][int(item_distance)]
                features[pdb_item].append([pdb_ppi_probability, ppi_probability, min_distance_frequency, distance_delta])
                '''
        return pdb_ppi_result,features

    @staticmethod
    def __rescore(pdb_ppi_result,features,data_package):
        score_dictionary={}
        for pdb_item in pdb_ppi_result:
            score_dictionary[pdb_item]=[]
            for order,item in enumerate(features[pdb_item]):
                if item:
                    item_ppi=pdb_ppi_result[pdb_item][order][0]
                    raw_score=min(data_package.my_CxResult.CFLOW6_SCORE[data_package.my_CxResult.CFLOW6_PPI.index(item_ppi)])
                    if item[0]*item[1]*(item[2]/(item[2]+1))*item[3]==0:
                        new_score=raw_score
                    else:
                        new_score=raw_score/(item[0]*item[1]*(item[2]/(item[2]+1))*item[3])
                    score_dictionary[pdb_item].append(new_score)
                else:
                    score_dictionary[pdb_item].append([])
        return score_dictionary

    @staticmethod
    def __rankAndOutPut(pdb_ppi_result,score_dictionary,data_package):
        '''
        score_list,pdb_ppi_pointer=[],[]
        for pdb_item in score_dictionary:
            for order,dis in enumerate(score_dictionary[pdb_item]):
                if dis==[]:
                    continue
                else:
                    score_list.append(dis)
                    pdb_ppi_pointer.append((pdb_item,order))
        score_list=(np.array(score_list)-np.min(score_list))/(np.max(score_list)-np.min(score_list))
        norm_score_list=10*(1-np.sin((np.pi/(2*(max(score_list)+0.1)))*np.array(score_list)))**2
        outputfilter_path = os.path.join(data_package.my_config.E1_PATH_EXPORT, f'structure_alignment_result_rerank_filtered.{data_package.my_config.C66_OUTPUT_DATE}.csv')
        output_path = os.path.join(data_package.my_config.E1_PATH_EXPORT, f'structure_alignment_result_rerank.{data_package.my_config.C66_OUTPUT_DATE}.csv')
        f2=open(output_path,'w')
        with open(outputfilter_path, 'w') as f:
            f.write('PDB code,PPI,Model,Protein 1 Chain,Protein 2 Chain,CA Distance,Score,Dist Command\n')
            f2.write('PDB code,PPI,Model,Protein 1 Chain,Protein 2 Chain,CA Distance,Score,Dist Command\n')
            for order in np.argsort(norm_score_list)[::-1]:
                if norm_score_list[order]<=9:
                    continue
                pdb_code,order_num=pdb_ppi_pointer[order]
                info_item=pdb_ppi_result[pdb_code][order_num]
                f2.write(
                    '{},{},{},{},{},{:.3f},{:.3f},"{}"\n'.format(pdb_code, info_item[0], info_item[2], info_item[3],
                                                                  info_item[4], info_item[5], 10(norm_score_list[order]-9),
                                                                  info_item[6]))
                if info_item[5] >= data_package.my_config.C65_LINKER_MAX_ARM_LENGTH:
                    continue
                f.write('{},{},{},{},{},{:.3f},{:.3f},"{}"\n'.format(pdb_code,info_item[0],info_item[2],info_item[3],info_item[4],info_item[5],10(norm_score_list[order]-9),info_item[6]))
        f2.close()
        '''
        rankAndOutput(pdb_ppi_result,score_dictionary,data_package)