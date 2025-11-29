import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)
sns.set_style("dark")
sns.set_style("ticks")
colorlist = ["#006662", "#006694", "#7a5195", "#ef5675", "#ffa600"]  # light

# MEB
algo1 = ["MBCI", "MBC", "BCviz", "BCviz+", "BCviz-"]
algo_lab1 = ["$MBCI$", "$MBC$", "$BCviz$", "$BCviz+$", "$BCviz-$"]
# MVB
algo2 = ["BIIP", "MBC", "BCviz", "BCviz+", "BCviz-"]
algo_lab2 = ["$BIIP$", "$MBC$", "$BCviz$", "$BCviz+$", "$BCviz-$"]
# MBB
algo3 = ["ooMBC", "MBC", "BCviz", "BCviz+", "BCviz-"]
algo_lab3 = ["$ooMBB$", "$MBC$", "$BCviz$", "$BCviz+$", "$BCviz-$"]

datasets = ["marvel", "writer", "BookCrossing",
            "Team", "Actor-movie", "Twitter",
            "WikiPedia", "DBLP"]
datasets_lab = ["$MV$", "$WT$", "$BC$",
            "$TM$", "$AM$", "$TT$",
            "$WP$", "$DB$"]


def read_data(path, pro_type, algo_name):
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
        if datalist[0] == pro_type and datalist[1] == algo_name:
            maxy = max(maxy, float(datalist[5]))  # search
            y.append(float(datalist[5]))  # ms
            x.append(row)
            # Csize.append(int(datalist[6]))
            row = row + 1
    return x, y, Csize, maxy


def Draw5(rcount, x_label, title, labels):
    plt.figure(figsize=(25, 10))
    plt.rcParams['figure.figsize'] = (80, 40)
    plt.xlabel("", fontsize=60)  # 49 50
    plt.xticks(fontsize=48)  # 28 40
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=True
    )
    plt.ylabel("Search Time(ms)", fontsize=60)  # 49 80
    plt.yticks(fontsize=60)  # 25 52 40
    plt.tick_params(axis='y', width=2, length=20)
    plt.title(title, size=60)
    rwidth = 2.5

    x1 = [0, 15, 30, 45, 60, 75, 90, 105]
    x2 = [i + rwidth for i in x1]
    x0 = [i + rwidth*1.5 for i in x1]
    x3 = [i + rwidth*2 for i in x1]
    x4 = [i + rwidth*3 for i in x1]
    x5 = [i + rwidth*4 for i in x1]

    y1 = [i[0] for i in rcount]
    y2 = [i[1] for i in rcount]
    y3 = [i[2] for i in rcount]
    y4 = [i[3] for i in rcount]
    y0 = [0, 0, 0, 0, 0, 0, 0, 0]
    y5 = [i[4] for i in rcount]

    b1 = plt.bar(x1, y1, width=rwidth, color=colorlist[0], alpha=0.8, label=labels[0], log=1)
    b2 = plt.bar(x2, y2, width=rwidth, color=colorlist[1], alpha=0.8, label=labels[1], log=1)
    b3 = plt.bar(x3, y3, width=rwidth, color=colorlist[2], alpha=0.8, label=labels[2], log=1, tick_label=x_label)
    b4 = plt.bar(x4, y4, width=rwidth, color=colorlist[3], alpha=0.8, label=labels[3], log=1)
    b5 = plt.bar(x5, y5, width=rwidth, color=colorlist[4], alpha=0.8, label=labels[4], log=1)

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

    plt.legend(handles=[b1, b2, b3, b4, b5], loc='best', fontsize=45, ncol=3)  # 23 40 30 frameon=False,

    # plt.xticks(rotation=-15)  # -18
    plt.xticks()
    # plt.legend()
    file = title + "-search-time.jpg"
    plt.savefig(file)
    print(file + " has been generated.")
    # plt.show()


if __name__ == '__main__':
    print("---------Search Time----------")
    path_algo = "search-time.txt"
    path_baseline = "baseline-search-time.txt"
    problems = ["MEB", "MVB", "MBB"]
    algo = [algo1, algo2, algo3]
    algo_lab = [algo_lab1, algo_lab2, algo_lab3]

    for i in range(0, len(problems)):
        pro_type = problems[i]
        print("Problem type:", pro_type)
        arrayY = []
        maxy = 0

        for j in range(0, len(algo[i])):
            al = algo[i][j]
            path = ""
            if j <= 1:
                path = path_baseline
            else:
                path = path_algo
            # print(al, path)
            X, Y, Csize_array, up = read_data(path, pro_type, al)
            arrayY.append(Y)
            maxy = max(maxy, up)
        array = np.transpose(arrayY)  # array.T
        Draw5(array, datasets_lab, pro_type, algo_lab[i])

        # timesBBVS = 0.0
        # timesSBVS = 0.0
        # timesBCviz = 0.0
        # for i in range(0, len(array[0])):
        #     timesBCviz = max(array[i][0] / array[i][1], timesSBVS)
        #     timesBBVS = max(array[i][0] / array[i][2], timesBBVS)
        #     timesSBVS = max(array[i][0] / array[i][3], timesSBVS)
        # print(timesBCviz)
        # print(timesBBVS)
        # print(timesSBVS)
        # print(timesSBVS / timesBBVS)

    # tall 0.6 right 0.9
    # 2308x877
