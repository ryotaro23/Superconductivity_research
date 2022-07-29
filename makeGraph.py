import pandas as pd
import matplotlib.pyplot as plt
import os
import ntpath
import traceback
import glob
import pathlib

def get_cycle_info(df):
    try:
        cycle_start = min(list(df['a']))
        cycle_end = max(list(df['a']))
        dx = round(df['a'].iloc[1] - df['a'].iloc[0],4)
        one_cycle = int((cycle_end - cycle_start)/dx) + 1
        Plus = max(list(df['b']))
        Minus = min(list(df['b']))
        abs_Plus = abs(Plus)
        abs_Minus = abs(Minus)
        if abs_Plus < abs_Minus:
            Abs = Minus
        elif abs_Plus > abs_Minus:
            Abs = Plus

    except Exception as e:
        print(traceback.format_exc())

    return {'start':cycle_start,'end':cycle_end,'one_cycle':one_cycle,'dx':dx,'abs_max':Abs}

def main():
    try:
        #番号を取得
        num_path = "./sim1/num_data_save.txt"
        num = int(open(num_path,'r', encoding='UTF-8').read())-1

        dir_path = "./sim1/VRE_MD_140/VRE_MD_file_"
        folder_name = "/VRE_MD_140_"
        file_path = dir_path + str(num) + folder_name + str(num)
        save_fig_path = dir_path + str(num) + "/figures_VRE_MD_140_" + str(num) + "/"

        #グラフを保存するフォルダーを作成
        p = pathlib.Path(save_fig_path)
        p.mkdir(exist_ok=True)

        file = file_path + ".txt"

        data = open(file, 'r').readlines()
        data_list = []

        #txtファイルの中身を全てデータフレームに格納
        for d in data:
            data_list.append(d.split(','))

        #データフレームに格納された数値を全て小数にする
        df = pd.DataFrame(data_list,columns = ['a','b','c','d','e','Rl','g','Rs','Dp_L2M','j','Dp_S2L']) #右から1,3,4,6
        df['a'] = df['a'].astype(float)
        df['b'] = df['b'].astype(float)
        df['c'] = df['c'].astype(float)
        df['d'] = df['d'].astype(float)
        df['e'] = df['e'].astype(float)
        df['Rl'] = df['Rl'].astype(float)
        df['g'] = df['g'].astype(float)
        df['Rs'] = df['Rs'].astype(float)
        df['Dp_L2M'] = df['Dp_L2M'].astype(float)
        df['j'] = df['j'].astype(float)
        df['Dp_S2L'] = df['Dp_S2L'].str.replace('\n','')
        df['Dp_S2L'] = df['Dp_S2L'].astype(float)

        #サイクルの情報
        one_cycle = get_cycle_info(df)['one_cycle']
        total_cycle = int(len(data)/one_cycle)

        #deta_listを空にする
        data_list.clear()

        #サイクルの数だけグラフを作る
        #平均のグラフの絶対値が一番大きいところが出力されるようにする
        for n in range(total_cycle):
            df1 = df[:one_cycle]
            df2 = df[one_cycle:]
            abs_max = get_cycle_info(df1)['abs_max']
            Rl = df1['Rl'].iloc[0]
            Rs = df1['Rs'].iloc[0]
            L2M = df1['Dp_L2M'].iloc[0]
            S2L = df1['Dp_S2L'].iloc[0]
            fig = plt.figure()
            plt.scatter(df1['a'],df1['b'])
            plt.scatter(df1['a'],df1['c'])
            plt.scatter(df1['a'],df1['d'])
            graph_title = ntpath.basename(file).replace('.txt','')
            plt.title('{0}(b={1},Rl={2},Rs={3},Dp_L2M={4},Dp_S2L={5})'.format(graph_title,abs_max,Rl,Rs,L2M,S2L)) #Rl=6番目、Rs=4番目、Dp_L2M=3番目、Dp_S2L = 1番目
            plt.xlabel('Fl/Fp')
            plt.ylabel('V_Ave')
            plt.show()
            fig.savefig('{0}(b={1},Rl={2},Rs={3},Dp_L2M={4},Dp_S2L={5}).png'.format(save_fig_path,abs_max,Rl,Rs,L2M,S2L))
            df = df.drop(range(one_cycle))
            df = df.reset_index(drop=True)

    #エラーはエラーコードを表示してスルーする
    except Exception as e:
        print(traceback.format_exc())
        return

#実行
main()
