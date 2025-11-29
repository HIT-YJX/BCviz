import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)
sns.set_style("dark")
sns.set_style("ticks")
path = "./"

# colorlist = ["#093034", "#7a5195", "#ef5675", "#ffa600"]  # best
colorlist = ["#7a5195", "#ef5675", "#ffa600"]  # light
colorlist4 = ["#006662", "#7a5195", "#ef5675", "#ffa600"]  # light

algo = ["BCviz", "BCviz+", "BCviz-"]
algo4 = ["MBCI", "BCviz", "BCviz+", "BCviz-"]
algo_lab = ["$BCviz$", "$BCviz+$", "$BCviz-$"]
algo_lab4 = ["$MBCI$", "$BCviz$", "$BCviz+$", "$BCviz-$"]
datasets = ["marvel", "writer", "BookCrossing",
            "Team", "Actor-movie", "Twitter",
            "WikiPedia", "DBLP"]
datasets_lab = ["$MV$", "$WT$", "$BC$",
            "$TM$", "$AM$", "$TT$",
            "$WP$", "$DB$"]


def read_data(path, sizetype, algo_name):
    file = open(path, "r")
    x = []
    y = []
    Csize = []
    row = 0
    maxy = 0
    for i, d in enumerate(file):
        datalist = d.strip().split()
        if len(datalist) < 6:
            continue
        if datalist[0] == sizetype:
            if datalist[1] == algo_name:
                maxy = max(maxy, float(datalist[5]))  # index time
                y.append(float(datalist[5]))
                x.append(row)
                # Csize.append(int(datalist[7]))
                row = row + 1
    return x, y, Csize, maxy


def Draw3(rcount, x_label, title):
    plt.figure(figsize=(25, 10))
    plt.rcParams['figure.figsize'] = (80, 50)
    plt.xlabel("", fontsize=60)
    plt.xticks(fontsize=48)  # 28
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=True
    )
    plt.ylabel("Construction Time(s)", fontsize=60)
    plt.yticks(fontsize=60)  # 25
    plt.tick_params(axis='y', width=2, length=18)

    rwidth = 3.7  # 4 3.7 3
    x1 = [0, 15, 30, 45, 60, 75, 90, 105]
    x2 = [i + rwidth*1 for i in x1]
    x0 = [i + rwidth*1.5 for i in x1]
    x3 = [i + rwidth*2 for i in x1]

    y1 = [i[0] for i in rcount]
    y2 = [i[1] for i in rcount]
    y3 = [i[2] for i in rcount]
    y0 = [0, 0, 0, 0, 0, 0, 0, 0]

    alpha1 = 0.8
    b1 = plt.bar(x1, y1, width=rwidth, color=colorlist[0], alpha=alpha1, label=algo_lab[0], log=1)
    b2 = plt.bar(x2, y2, width=rwidth, color=colorlist[1], alpha=alpha1, label=algo_lab[1], log=1, tick_label=x_label)
    b3 = plt.bar(x3, y3, width=rwidth, color=colorlist[2], alpha=alpha1, label=algo_lab[2], log=1)

    ax = plt.gca()
    # set border color
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    # set border width-size
    bwith = 2  # 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    plt.legend(handles=[b1, b2, b3], loc='best', fontsize=45, ncol=1)  # 35 23 40 30 frameon=False,
    plt.title(title, size=60)
    filename = title + "-IndexTime.jpg"
    print(filename + " has been generated.")
    plt.savefig(filename)
    # plt.show()


def Draw4(rcount, x_label, title):
    plt.figure(figsize=(25, 10))
    plt.rcParams['figure.figsize'] = (80, 50)
    plt.xlabel("", fontsize=60)
    plt.xticks(fontsize=48)  # 28
    plt.tick_params(
        axis='x',  # Specify the axis to adjust, 'x' or 'y'.
        which='both',  # Adjust the position of tick labels.
        bottom=False,  # Control the display of bottom tick labels.
        top=False,  # Control the display of top tick labels.
        labelbottom=True  # Control the text display of bottom tick labels.
    )
    plt.ylabel("Construction Time(s)", fontsize=60)
    plt.yticks(fontsize=60)  # 25
    plt.tick_params(axis='y', width=2, length=18)

    rwidth = 3  # 4 3.7 3
    x1 = [0, 15, 30, 45, 60, 75, 90, 105]
    x2 = [i + rwidth*1 for i in x1]
    x0 = [i + rwidth*1.5 for i in x1]
    x3 = [i + rwidth*2 for i in x1]
    x4 = [i + rwidth*3 for i in x1]

    y1 = [i[0] for i in rcount]
    y2 = [i[1] for i in rcount]
    y3 = [i[2] for i in rcount]
    y4 = [i[3] for i in rcount]
    y0 = [0, 0, 0, 0, 0, 0, 0, 0]

    alpha1 = 0.8
    b1 = plt.bar(x1, y1, width=rwidth, color=colorlist4[0], alpha=alpha1, label=algo_lab4[0], log=1)
    b2 = plt.bar(x2, y2, width=rwidth, color=colorlist4[1], alpha=alpha1, label=algo_lab4[1], log=1)  # , tick_label=x_label
    b0 = plt.bar(x0, y0, width=0.1, color='white', tick_label=x_label)
    b3 = plt.bar(x3, y3, width=rwidth, color=colorlist4[2], alpha=alpha1, label=algo_lab4[2], log=1)
    b4 = plt.bar(x4, y4, width=rwidth, color=colorlist4[3], alpha=alpha1, label=algo_lab4[3], log=1)

    ax = plt.gca()
    # set border color
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    # set border width-size
    bwith = 2  # 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    plt.legend(handles=[b1, b2, b3, b4], loc='best', fontsize=45, ncol=2)  # 35 23 40 30 frameon=False,
    plt.title(title, size=60)
    filename = title + "-IndexTime.jpg"
    plt.savefig(filename)
    print(filename + " has been generated.")
    # plt.show()


if __name__ == '__main__':
    print("---------Index Construction Time----------")
    path_algo = "index-time.txt"
    path_baseline = "baseline-index-time.txt"
    problems = ["MEB", "MVB", "MBB"]

    for pro in problems:
        print("Problem type:", pro)
        arrayY = []
        maxy = 0
        if pro == "MEB":
            X, Y, Csize_array, up = read_data(path_baseline, pro, "MBCI")
            arrayY.append(Y)
            maxy = max(maxy, up)
        for al in algo:
            X, Y, Csize_array, up = read_data(path_algo, pro, al)
            arrayY.append(Y)
            maxy = max(maxy, up)
        array = np.transpose(arrayY)  # array.T
        rateBCviz1toBCviz = []
        rateBCviz2toBCviz = []
        if pro == "MEB":
            Draw4(array, datasets_lab, pro)
            for i in range(0, len(array)):
                rateBCviz1toBCviz.append(array[i][1] / array[i][2])
                rateBCviz2toBCviz.append(array[i][2] / array[i][3])
            print("Table 4. Speedup of index construction (MEB, s_min = t_min = 3).")
            print("BCviz/BCviz+:", rateBCviz1toBCviz)
            print("BCviz+/BCviz-:", rateBCviz2toBCviz)
        else:
            Draw3(array, datasets_lab, pro)
    # tall size 0.81 & 0.76 & tall 0.6 right 0.9
    # tall 1092 & 1024 & 2308x877
