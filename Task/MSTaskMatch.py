import os

from System.MSData import CDataPack
from Function.MSFunctionIO import CFunctionParsePDB, CFunctionParseCIF, CFunctionOutPut
from Function.MSFunctionStructure import CFunctionSequenceAlignment, CFunctionPPiPDB
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from itertools import repeat


class CTaskSequenceAlignment:

    def work_pdb(self, data_package=CDataPack()):
        extract_function = CFunctionParsePDB()
        alignment_task = CFunctionSequenceAlignment()
        for pdb_name in data_package.my_config.C63_SUCCESS_PDBS:
            pdb_sequences,pdb_sequence_order = extract_function.parsePDBSequence(pdb_name,data_package.my_config.C60_TYPE01_DATA_PATH)
            alignment_dictionary = alignment_task.sequenceAlignment(pdb_sequences,pdb_sequence_order,
                                                                    data_package.my_config.C63_PDB_PROTEINS_MAP[
                                                                        pdb_name], data_package)
            data_package.my_config.C64_SEQUENCE_ALIGNMENT_RESULT[pdb_name] = alignment_dictionary
        os.remove(data_package.my_config.C60_TYPE01_DATA_PATH+r'\query_sequences.fasta')
        os.remove(data_package.my_config.C60_TYPE01_DATA_PATH+r'\subject_sequence.fasta')
        os.remove(data_package.my_config.C60_TYPE01_DATA_PATH+r'\result.xml')

    def work_cif(self, data_package=CDataPack()):
        index_list=list(range(len(data_package.my_config.C63_SUCCESS_PDBS)))
        with ProcessPoolExecutor(max_workers=data_package.my_config.P1_NUMBER_THREAD) as pool:
            results = list(pool.map(self.single_cif_worker,data_package.my_config.C63_SUCCESS_PDBS,
                                    repeat(data_package.my_config.C63_PDB_PROTEINS_MAP),
                                    repeat(data_package.my_ProteinIndex.CFLOW6_SEQUENCE_INDEX),
                                    repeat(data_package.my_config.C60_TYPE01_DATA_PATH),
                                    repeat(data_package.my_config.C64_SEQUENCE_ALIGNMENT_E_VALUE),
                                    index_list))

        for i in index_list:
            os.remove(data_package.my_config.C60_TYPE01_DATA_PATH + fr'\query_sequences.{i}.fasta')
            os.remove(data_package.my_config.C60_TYPE01_DATA_PATH + fr'\subject_sequence.{i}.fasta')
            os.remove(data_package.my_config.C60_TYPE01_DATA_PATH + fr'\result.{i}.xml')
        for pdb_name,result_item in zip(data_package.my_config.C63_SUCCESS_PDBS,results):
            data_package.my_config.C64_SEQUENCE_ALIGNMENT_RESULT[pdb_name] = result_item


    def single_cif_worker(self,pdb_name,dp_protein_map,dp_seq_index,dp_data_path,dp_e_value,index):
        extract_function = CFunctionParseCIF()
        alignment_task = CFunctionSequenceAlignment()
        #pdb_sequences, pdb_sequence_order = extract_function.parseCIFSequence(pdb_name,dp_data_path)
        pdb_sequences, pdb_sequence_order = extract_function.parseCIFSequenceFromAtom(pdb_name, dp_data_path)
        if pdb_sequences:
            alignment_dictionary = alignment_task.sequenceAlignment(pdb_sequences, pdb_sequence_order,
                                                                dp_protein_map[pdb_name], dp_seq_index,dp_data_path,dp_e_value,index)
        else:
            alignment_dictionary={}
        return alignment_dictionary

class CTaskStructureDistance:

    def work_pdb(self, data_package=CDataPack()):
        extract_function = CFunctionParsePDB()
        distance_function = CFunctionPPiPDB()
        for ppi_index in range(len(data_package.my_CxResult.CFLOW6_PPI)):
            data_package.my_CxResult.CFLOW6_DISTANCE.append({})
            for pdb_structure in data_package.my_CxResult.CFLOW6_PPI_PDB[ppi_index]:
                if pdb_structure not in data_package.my_config.C63_SUCCESS_PDBS:
                    continue
                ppi_item = data_package.my_CxResult.CFLOW6_PPI[ppi_index]
                coordinates = extract_function.parsePDBCoordinate(pdb_structure,data_package.my_config.C60_TYPE01_DATA_PATH)
                return_distance=distance_function.calculatePPiPDBDistance(pdb_structure, ppi_item.protein1, ppi_item.protein2,
                                                          ppi_item.site1, ppi_item.site2, coordinates, data_package)
                data_package.my_CxResult.CFLOW6_DISTANCE[-1][pdb_structure]=return_distance

    def work_cif(self, data_package=CDataPack()):
        task_list=[]
        for index, ppi in enumerate(data_package.my_CxResult.CFLOW6_PPI):
            for pdb in data_package.my_CxResult.CFLOW6_PPI_PDB[index]:
                task_list.append((index,ppi.protein1,ppi.protein2,ppi.site1,ppi.site2,pdb))
        for _ in range(len(data_package.my_CxResult.CFLOW6_PPI)):
            data_package.my_CxResult.CFLOW6_DISTANCE.append({})
        each_processor_task_number=len(task_list)//data_package.my_config.P1_NUMBER_THREAD
        range_box=[[i*each_processor_task_number,(i+1)*each_processor_task_number] for i in range(data_package.my_config.P1_NUMBER_THREAD)]
        range_box[-1][1]=len(task_list)
        submit_task_list=[task_list[item[0]:item[1]] for item in range_box]
        with ProcessPoolExecutor(max_workers=data_package.my_config.P1_NUMBER_THREAD) as pool:
            results = list(pool.map(self.single_cif_worker,submit_task_list,
                                    repeat(data_package.my_config.C64_SEQUENCE_ALIGNMENT_RESULT),
                                    repeat(data_package.my_ProteinIndex.CFLOW6_SEQUENCE_INDEX),
                                    repeat(data_package.my_config.C63_SUCCESS_PDBS),
                                    repeat(data_package.my_config.C60_TYPE01_DATA_PATH)))
        flatten_results=[item2 for item in results for item2 in item]
        for item in flatten_results:
            if item:
                idx,result_list=item
                if result_list:
                    data_package.my_CxResult.CFLOW6_DISTANCE[idx][result_list[0][0]] = result_list


    def single_cif_worker(self,submit_task_list,alignment_result,sequence_index,success_pdbs,data_path):
        return_list=[]
        for task_item in submit_task_list:
            index,ppi_p1, ppi_p2, ppi_s1, ppi_s2, pdb_structure=task_item
            extract_function = CFunctionParseCIF()
            distance_function = CFunctionPPiPDB()
            if pdb_structure not in success_pdbs:
                return_list.append((index,[]))
                continue
            #coordinates = extract_function.parseCIFCoordinate(pdb_structure,data_path)
            coordinates = extract_function.parseCIFCoordinateFromAtom(pdb_structure, data_path)
            if coordinates:
                return_distance = distance_function.calculatePPiPDBDistance(pdb_structure,ppi_p1,
                                                                            ppi_p2,ppi_s1, ppi_s2, coordinates,
                                                                            alignment_result,sequence_index)
            else:
                return_distance=[]
            return_list.append((index,return_distance))
        return return_list