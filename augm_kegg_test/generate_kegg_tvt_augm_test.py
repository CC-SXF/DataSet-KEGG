# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 19:47:08 2022

@author: CC-SXF
"""

import re
# import sys
import numpy as np
from rdkit import Chem
# from itertools import permutations
from os import listdir, mkdir, path



if 'datas_sfcv_augm' not in listdir('datas'):
    mkdir('./datas/datas_sfcv_augm/')

dir_sour = './datas/datas_sfcv/'
dir_augm = './datas/datas_sfcv_augm/'

# create folders
subdatas = ['SFCV_1_1', 'SFCV_2_1', 'SFCV_3_1', 'SFCV_4_1', 'SFCV_5_1']
for subdata in subdatas:
    if subdata not in listdir(dir_augm):
        mkdir(path.join(dir_augm, subdata))


def smi_tokenizer(smi):
    """
        # Tokenize a SMILES molecule or reaction
        # https://github.com/pschwllr/MolecularTransformer
    """
    # import re
    pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)


def enumsmi(smi, n_augm=10, postfix=None, random_seed=0):
    """
         SMILES enumeration
    """
    np.random.seed(random_seed)
    mol = Chem.MolFromSmiles(smi)
    n_at = mol.GetNumAtoms()
    m_neworder = list(range(n_at))

    smi_temp = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    assert smi_temp == smi

    augm_smi_list = [smi_temp]
    n_augm_temp = 1

    for n in range(n_augm*100):
        if n_augm_temp >= n_augm:
            break
        #
        nm = Chem.RenumberAtoms(mol, m_neworder)
        smi_temp = Chem.MolToSmiles(nm, canonical=False)
        if smi_temp not in augm_smi_list:
            augm_smi_list.append(smi_temp)
            n_augm_temp += 1
        #
        np.random.shuffle(m_neworder)

    # 分词
    augm_smi_list = [smi_tokenizer(smi_temp) for smi_temp in augm_smi_list]

    # 补充后缀
    if postfix != None:
        augm_smi_list = [' > '.join([smi_temp, postfix]) for smi_temp in augm_smi_list]

    # 随机补全
    num_random_choice = (n_augm - len(augm_smi_list))
    if num_random_choice != 0:
        augm_smi_list += list(np.random.choice(augm_smi_list, num_random_choice))

    return augm_smi_list



def augmdata(subdata='SFCV_1_1', n_augm=10, dir_augm = './datas/datas_sfcv_augm/', dir_sour='./datas/datas_sfcv/'):
    """
    """
    # test
    with open(path.join(dir_sour, subdata, 'src-test.txt'), 'r') as test_sour_src_file:
        test_sour_src_lines = test_sour_src_file.readlines()
    with open(path.join(dir_sour, subdata, 'tgt-test.txt'), 'r') as test_sour_tgt_file:
        test_sour_tgt_lines = test_sour_tgt_file.readlines()
    test_sour_src_lines = [item.strip() for item in test_sour_src_lines]
    test_sour_tgt_lines = [item.strip() for item in test_sour_tgt_lines]

    # print(len(test_sour_src_lines))
    # print(len(test_sour_tgt_lines))
    print('source:', test_sour_src_lines[0])
    print('target:', test_sour_tgt_lines[0])
    print('')

    # test source (smiles, enzyme)
    test_sour_src_smi_list = []
    test_sour_src_ezy_list = []
    for line in test_sour_src_lines:
        (smi, ezy) = line.split(' > ')
        smi = smi.replace(' ', '')
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        assert smi != ''
        test_sour_src_smi_list.append(smi)
        test_sour_src_ezy_list.append(ezy)

    # print(len(test_sour_src_smi_list))
    # print(len(test_sour_src_ezy_list))
    print('source:', test_sour_src_smi_list[0])
    # print(test_sour_src_ezy_list[0])

    # test target (smiles)
    test_sour_tgt_smi_list = []
    for line in test_sour_tgt_lines:
        smi = line.replace(' ', '')
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        assert smi != ''
        test_sour_tgt_smi_list.append(smi)

    # print(len(test_sour_tgt_smi_list))
    print('target:', test_sour_tgt_smi_list[0])


    # augmentation: augm test source (smiles, enzyme)
    test_sour_src_augmsmi_list = []
    for smi, ezy in zip(test_sour_src_smi_list, test_sour_src_ezy_list):
        augmsmi = enumsmi(smi, n_augm, ezy)
        test_sour_src_augmsmi_list.append(augmsmi)

    # augmentation: test target (smiles)
    test_sour_tgt_augmsmi_list = []
    for smi in test_sour_tgt_smi_list:
        """
        # # augmsmi = enumsmi(smi, n_augm)
        """
        # ########################################################################################
        augmsmi = [smi_tokenizer(smi)] * n_augm
        # ########################################################################################
        test_sour_tgt_augmsmi_list.append(augmsmi)

    """
    # 初始化
    augm_src_test = [[] for idx in range(n_augm)]
    aubm_tgt_test = [[] for idx in range(n_augm)]

    # 数据重构
    for src_augmsmi, tgt_augmsmi in zip(test_sour_src_augmsmi_list, test_sour_tgt_augmsmi_list):
        for idx in range(n_augm):
            augm_src_test[idx].append(src_augmsmi[idx])
            aubm_tgt_test[idx].append(tgt_augmsmi[idx])
    # return augm_src_test, aubm_tgt_test
    """

    augm_src_test = test_sour_src_augmsmi_list
    aubm_tgt_test = test_sour_tgt_augmsmi_list


    # ‘\n’ 训练子数据, 内部连接
    augm_src_test = ['\n'.join(item) for item in augm_src_test]
    aubm_tgt_test = ['\n'.join(item) for item in aubm_tgt_test]
    # return augm_src_test, aubm_tgt_test

    # ‘\n’ 训练子数据, 之间连接
    augm_src_test_content = '\n'.join(augm_src_test)
    aubm_tgt_test_content = '\n'.join(aubm_tgt_test)
    # return augm_src_test_content, aubm_tgt_test_content


    # test (增强之后的数据写入文件)
    with open(path.join(dir_augm, subdata, 'src-test.txt'), 'w') as test_augm_src_file:
        test_augm_src_file.write(augm_src_test_content)
    with open(path.join(dir_augm, subdata, 'tgt-test.txt'), 'w') as test_augm_tgt_file:
        test_augm_tgt_file.write(aubm_tgt_test_content)

    # valid
    with open(path.join(dir_sour, subdata, 'src-valid.txt'), 'r') as valid_sour_src_file:
        content = valid_sour_src_file.read()
        with open(path.join(dir_augm, subdata, 'src-valid.txt'), 'w') as valid_augm_src_file:
            valid_augm_src_file.write(content)
    with open(path.join(dir_sour, subdata, 'tgt-valid.txt'), 'r') as valid_sour_tgt_file:
        content = valid_sour_tgt_file.read()
        with open(path.join(dir_augm, subdata, 'tgt-valid.txt'), 'w') as valid_augm_tgt_file:
            valid_augm_tgt_file.write(content)
    # train
    with open(path.join(dir_sour, subdata, 'src-train.txt'), 'r') as train_sour_src_file:
        content = train_sour_src_file.read()
        with open(path.join(dir_augm, subdata, 'src-train.txt'), 'w') as train_augm_src_file:
            train_augm_src_file.write(content)
    with open(path.join(dir_sour, subdata, 'tgt-train.txt'), 'r') as train_sour_tgt_file:
        content = train_sour_tgt_file.read()
        with open(path.join(dir_augm, subdata, 'tgt-train.txt'), 'w') as train_augm_tgt_file:
            train_augm_tgt_file.write(content)
    #
    pass



def augm_valid(subdata='SFCV_1_1', n_augm=10, dir_augm = './datas/datas_sfcv_augm/', dir_sour='./datas/datas_sfcv/'):
    """
    """

    with open(path.join(dir_sour, subdata, 'src-test.txt'), 'r') as test_sour_src_file:
        test_sour_src_lines = test_sour_src_file.readlines()
    with open(path.join(dir_sour, subdata, 'tgt-test.txt'), 'r') as test_sour_tgt_file:
        test_sour_tgt_lines = test_sour_tgt_file.readlines()
    test_sour_src_lines = [item.strip() for item in test_sour_src_lines]
    test_sour_tgt_lines = [item.strip() for item in test_sour_tgt_lines]
    #
    # test source (smiles, enzyme)
    test_sour_src_smi_list = []
    for line in test_sour_src_lines:
        (smi, _) = line.split(' > ')
        smi = smi.replace(' ', '')
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        assert smi != ''
        test_sour_src_smi_list.append(smi)
    #
    # test target (smiles)
    test_sour_tgt_smi_list = []
    for line in test_sour_tgt_lines:
        smi = line.replace(' ', '')
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        assert smi != ''
        test_sour_tgt_smi_list.append(smi)


    with open(path.join(dir_augm, subdata, 'src-test.txt'), 'r') as test_augm_src_file:
        test_augm_src_lines = test_augm_src_file.readlines()
    with open(path.join(dir_augm, subdata, 'tgt-test.txt'), 'r') as test_augm_tgt_file:
        test_augm_tgt_lines = test_augm_tgt_file.readlines()
    test_augm_src_lines = [item.strip() for item in test_augm_src_lines]
    test_augm_tgt_lines = [item.strip() for item in test_augm_tgt_lines]
    #
    # test source (smiles, enzyme)
    test_augm_src_smi_list = []
    for line in test_augm_src_lines:
        (smi, _) = line.split(' > ')
        smi = smi.replace(' ', '')
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        assert smi != ''
        test_augm_src_smi_list.append(smi)
    #
    # test target (smiles)
    test_augm_tgt_smi_list = []
    for line in test_augm_tgt_lines:
        smi = line.replace(' ', '')
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        assert smi != ''
        test_augm_tgt_smi_list.append(smi)

    print('')
    print('xN:', n_augm)
    print('src-test:', len(test_sour_src_smi_list))  # 2541
    print('tgt-test:', len(test_sour_tgt_smi_list))  # 2541
    print('src-test xN:', len(test_augm_src_smi_list))
    print('tgt-test xN:', len(test_augm_tgt_smi_list))
    print('')

    """
    n_test = len(test_sour_src_lines)
    """

    for idx, smi in enumerate(test_augm_src_smi_list):
        """
        if smi != test_sour_src_smi_list[idx%n_test]:
            print('src:', idx%n_test)
            print(smi)
            print(test_sour_src_smi_list[idx%n_test])
            print("")
            # sys.exit()
        """
        if smi != test_sour_src_smi_list[idx//n_augm]:
            print('src:', idx//n_augm)
            print(smi)
            print(test_sour_src_smi_list[idx//n_augm])
            print("")

    for idx, smi in enumerate(test_augm_tgt_smi_list):
        """
        if smi != test_sour_tgt_smi_list[idx%n_test]:
            print('tgt:', idx%n_test)
            print(smi)
            print(test_sour_tgt_smi_list[idx%n_test])
            print("")
            # sys.exit()
        """
        if smi != test_sour_tgt_smi_list[idx//n_augm]:
            print('tgt:', idx//n_augm)
            print(smi)
            print(test_sour_tgt_smi_list[idx//n_augm])
            print("")
            # sys.exit()
    pass




if __name__ == '__main__':
    """
        # The pattern of atom symbol and atom-atom mapping number.
        # import re
        # pattern_AtSyMn = re.compile("\[\d*([a-zA-Z\*]+)@*[-\+]*\d*:(\d+)\]")

        x010, X050, x100
        cd /ext/PyTorch_Cc/augm_kegg/ && source activate && conda activate onmt12 && python generate_kegg_tvt_augm.py >> ./datas/datas_sfcv_augm.log 2>&1

    """
    # glucose:
    # enumsmi('OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O')
    # augmdata()

    n_augm = 2
    # n_augm = 5
    # n_augm = 10
    # n_augm = 50
    # n_augm = 100
    for subdata in subdatas:
        print(subdata, )
        augmdata(subdata, n_augm)
        augm_valid(subdata, n_augm)
        print('\n')
    pass





