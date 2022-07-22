import os
from xml.etree import ElementTree
from itertools import product
from math import sqrt
from random import SystemRandom
from System.MSData import CDataPack


class CFunctionPPiPDB:

    def calculatePPiPDB(self, data_package=CDataPack()):
        for ppi_item in data_package.my_CxResult.CFLOW6_PPI:
            if ppi_item.protein1 not in data_package.my_config.C63_PROTEINS_PDB_MAP or \
                    ppi_item.protein2 not in data_package.my_config.C63_PROTEINS_PDB_MAP:
                data_package.my_CxResult.CFLOW6_PPI_PDB.append([])
                continue
            protein1_PDB = set(data_package.my_config.C63_PROTEINS_PDB_MAP[ppi_item.protein1].index)
            protein2_PDB = set(data_package.my_config.C63_PROTEINS_PDB_MAP[ppi_item.protein2].index)
            data_package.my_CxResult.CFLOW6_PPI_PDB.append(list(protein1_PDB.intersection(protein2_PDB)))

    def calculatePPiPDBDistance(self, pdb_name, protein1, protein2, site1, site2, coordinates,
                                alignment_result,sequence_index):
        return_distance = []
        protein1_match_chains, protein2_match_chains = self.__getProteinMatchChain(pdb_name, protein1, protein2, site1,
                                                                                   site2, alignment_result,sequence_index)
        state = 0
        for protein1_state_chain, protein2_state_chain in zip(protein1_match_chains, protein2_match_chains):
            if protein1_state_chain and protein2_state_chain:
                for chains in product(protein1_state_chain, protein2_state_chain):
                    protein1_chain, protein2_chain = chains

                    protein1_coordinate = coordinates[state][protein1_chain][
                        alignment_result[pdb_name][(state, protein1_chain)][1][
                            site1 - 1]]
                    protein2_coordinate = coordinates[state][protein2_chain][
                        alignment_result[pdb_name][(state, protein2_chain)][1][
                            site2 - 1]]
                    distance = sqrt((float(protein1_coordinate[0]) - float(protein2_coordinate[0])) ** 2 +
                                    (float(protein1_coordinate[1]) - float(protein2_coordinate[1])) ** 2 +
                                    (float(protein1_coordinate[2]) - float(protein2_coordinate[2])) ** 2)
                    protein1_chain_name=alignment_result[pdb_name][(state, protein1_chain)][2][1]
                    protein2_chain_name=alignment_result[pdb_name][(state, protein2_chain)][2][1]
                    pdb_site1 = alignment_result[pdb_name][(state, protein1_chain)][1][site1 - 1]
                    pdb_site2 = alignment_result[pdb_name][(state, protein2_chain)][1][site2 - 1]
                    if len(protein1_coordinate)==4:
                        sub_name1=protein1_coordinate[3]
                    else:
                        sub_name1=''
                    if len(protein2_coordinate)==4:
                        sub_name2=protein2_coordinate[3]
                    else:
                        sub_name2=''
                    return_distance.append((pdb_name, state, protein1_chain_name, protein2_chain_name, distance,pdb_site1,pdb_site2,sub_name1,sub_name2))
            state += 1
        return return_distance

    def __getProteinMatchChain(self, pdb_name, protein1, protein2, site1, site2, alignment_result,sequence_index):
        protein1_match_chains = [[] for _ in range(1000)]
        protein2_match_chains = [[] for _ in range(1000)]
        for item in alignment_result[pdb_name]:
            if alignment_result[pdb_name][item][0] == protein1 and site1 - 1 in \
                    alignment_result[pdb_name][item][1] and \
                    alignment_result[pdb_name][item][1][site1 - 1] != 'Empty' and \
                    sequence_index[protein1][site1 - 1] == \
                    alignment_result[pdb_name][item][2][0][
                    alignment_result[pdb_name][item][3].index(alignment_result[pdb_name][item][1][site1 - 1])]:
                    protein1_match_chains[item[0]].append(item[1])
            if alignment_result[pdb_name][item][0] == protein2 and site2-1 in \
                    alignment_result[pdb_name][item][1] and \
                    alignment_result[pdb_name][item][1][site2-1] != 'Empty' and \
                    sequence_index[protein2][site2 - 1] == \
                    alignment_result[pdb_name][item][2][0][
                    alignment_result[pdb_name][item][3].index(alignment_result[pdb_name][item][1][site2 - 1])]:
                protein2_match_chains[item[0]].append(item[1])
        return protein1_match_chains, protein2_match_chains


class CFunctionSequenceAlignment:

    def sequenceAlignment(self, query_sequences,pdb_sequence_order, uniprot_names, dp_seq_index,dp_data_path,dp_e_value,index):
        for item in uniprot_names:
            if item not in dp_seq_index:
                print(item)
        query_seq=dp_data_path+fr'\query_sequences.{index}.fasta'
        subject_seq=dp_data_path+fr'\subject_sequence.{index}.fasta'
        result_f=dp_data_path+fr'\result.{index}.xml'
        blast_f=dp_data_path+fr'\blastp.exe'
        with open(query_seq, 'w') as query_temp_file:
            for order1, single_state in enumerate(query_sequences):
                for order2, single_query in enumerate(single_state):
                    query_temp_file.write(
                        '>lcl|query_sequence.{}.{}.|query sequence{}.{}\n'.format(order1, order2, order1, order2))
                    query_temp_file.write('{}\n'.format(single_query[0]))
        with open(subject_seq, 'w') as subject_temp_file:
            for order, uniprot_name in enumerate(uniprot_names):
                subject_temp_file.write('>lcl|{}|subject_sequence{}\n'.format(uniprot_name, order))
                subject_temp_file.write('{}\n'.format(dp_seq_index[uniprot_name]))
        cmd_str=f'{blast_f} -query {query_seq} -subject {subject_seq} -out {result_f} -outfmt 5 -max_hsps 1 -evalue {dp_e_value}'
        os.system(cmd_str)
        alignment_result_dictionary = {}
        with open(result_f, 'r') as f:
            tree = ElementTree.parse(f)
            for node in tree.iter():
                if node.tag == 'Iteration_query-def':
                    protein_match_state = int(node.text.split('|')[1].split('.')[1])
                    protein_match_chain = int(node.text.split('|')[1].split('.')[2])
                    match_protein_name = ''
                    continue
                if node.tag == 'Hit_id':
                    match_protein_name = node.text.split('|')[1]
                if node.tag == 'Hsp_qseq':
                    query_sequence = node.text
                if node.tag == 'Hsp_hseq':
                    subject_sequence = node.text
                if node.tag == 'Statistics':
                    if match_protein_name:
                        q_start_position = query_sequences[protein_match_state][protein_match_chain][0].index(
                            query_sequence.split('-')[0])
                        s_start_position = dp_seq_index[match_protein_name].index(
                            subject_sequence.split('-')[0])
                        alignment_result_dictionary[(protein_match_state, protein_match_chain)] = [match_protein_name,
                                                                                                   {}, query_sequences[
                                                                                                       protein_match_state][
                                                                                                       protein_match_chain],
                                                                                                   pdb_sequence_order[protein_match_state][protein_match_chain]]
                        for q, s in zip(query_sequence, subject_sequence):
                            if q == '-':
                                alignment_result_dictionary[(protein_match_state, protein_match_chain)][1][
                                    s_start_position] = 'Empty'
                                s_start_position += 1
                                continue
                            if s == '-':
                                q_start_position += 1
                                continue
                            alignment_result_dictionary[(protein_match_state, protein_match_chain)][1][
                                s_start_position] =  pdb_sequence_order[protein_match_state][protein_match_chain][q_start_position]
                            s_start_position += 1
                            q_start_position += 1

                    else:
                        alignment_result_dictionary[(protein_match_state, protein_match_chain)] = ['', {}, '',{}]
        return alignment_result_dictionary
