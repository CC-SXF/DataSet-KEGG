#!bin/sh

# screen -h 1000000 -S generate_kegg_tvt_augm_test          
# 2024
# cd /ext/PyTorch_Cc/augm_kegg_test/ && conda activate onmt12 && sh generate_kegg_tvt_augm_test.sh  && cd ~                   
# 2025
# cd /ext/PyTorch_Cc/augm_kegg_test/ && source activate && conda activate onmt12 && sh generate_kegg_tvt_augm_test.sh  && cd ~           

python generate_kegg_tvt_augm_test.py >> ./datas/datas_sfcv_augm_test.log 2>&1                            


# screen 语法
# screen [-AmRvx -ls -wipe][-d <作业名称>][-h <行数>][-r <作业名称>][-s <shell>][-S <作业名称>]
# 1.新建会话
# screen -h 100000 -S session_name 
# 2.列出所有会话
# screen -ls
# 3.断开会话
# screen -d session_name
# 4.重新连接会话（先断开再连接）
# screen -r session_name
# 5.清除dead会话（检查目前所有的screen作业，并删除已经无法使用的screen作业）
# screen -wipe
# 6.杀死已经断开的会话（谨慎使用！！！谨慎使用！！！谨慎使用！！！）
# screen -X -S session_name quit
# screen -X -S generate_kegg_tvt_augm_test quit                





