import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)
sns.set_style("dark")
sns.set_style("ticks")
colorlist = ["#006662", "#006694", "#7a5195", "#ef5675", "#ffa600"]
algo = ["MBCI", "MBC", "BCviz+", "BCviz-", "BCviz"]
algo_lab = ["$MBCI$", "$MBC$", "$BCviz+$", "$BCviz-$", "$BCviz$"]
tu = [3, 4, 5, 6]
tu_lab = ["$(3,3)$", "$(4,4)$", "$(5,5)$", "$(6,6)$"]


def read_data2(path, da, algo_name, pro):
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
        if datalist[0] == pro and datalist[2] == da and datalist[1] == algo_name:
            value = float(datalist[5])
            maxy = max(maxy, value)  # construct
            y.append(value)
        x.append(row)
        row = row + 1
    return x, y, Csize, maxy


def Draw5(rcount, x_label, title):
    plt.figure(figsize=(28, 17))
    plt.rcParams['figure.figsize'] = (100, 70)
    plt.xlabel("($s$, $t$)", fontsize=70)
    plt.xticks(fontsize=70)  # 28 48
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=True,
    )
    plt.ylabel("Search Time(ms)", fontsize=70)
    plt.yticks(fontsize=70)  # 25
    plt.tick_params(axis='y', width=4, length=30)

    rwidth = 5
    x1 = [0, 30, 60, 90]
    x2 = [i + rwidth*1 for i in x1]
    x3 = [i + rwidth*2 for i in x1]
    x0 = [i + rwidth*1.5 for i in x1]
    x4 = [i + rwidth*3 for i in x1]
    x5 = [i + rwidth*4 for i in x1]

    y1 = [i[0] for i in rcount]
    y2 = [i[1] for i in rcount]
    y3 = [i[2] for i in rcount]
    y4 = [i[3] for i in rcount]
    y0 = [0, 0, 0]
    y5 = [i[4] for i in rcount]

    alpha1 = 0.8
    b1 = plt.bar(x1, y1, width=rwidth, color=colorlist[0], alpha=alpha1, label=algo_lab[0], log=1)
    # x3[2] = x1[2] + rwidth * 2.5
    # x3[4] = x1[4] + rwidth * 1.5
    # x3[6] = x1[6] + rwidth * 1
    # b0 = plt.bar(x0, y0, width=0.1, color='white', tick_label=x_label)
    b2 = plt.bar(x2, y2, width=rwidth, color=colorlist[1], alpha=alpha1, label=algo_lab[1], log=1)
    b3 = plt.bar(x3, y3, width=rwidth, color=colorlist[2], alpha=alpha1, label=algo_lab[2], log=1, tick_label=x_label)
    b4 = plt.bar(x4, y4, width=rwidth, color=colorlist[3], alpha=alpha1, label=algo_lab[3], log=1)
    b5 = plt.bar(x5, y5, width=rwidth, color=colorlist[4], alpha=alpha1, label=algo_lab[4], log=1)

    ax = plt.gca()
    # set border color
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    # set border width-size
    bwith = 4  # 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    plt.title(title, size=70)
    plt.legend(handles=[b1, b2, b3, b4, b5], loc="upper left",  fontsize=55, ncol=1)  # frameon=False,
    file = title + "-varied-search-time.jpg"
    plt.savefig(file)
    print(file + " has been generated.")
    # plt.show()


if __name__ == '__main__':
    print("-----Search Time (varied parameters)-----")
    path_algo = "varied-search-time.txt"
    path_baseline = "baseline-varied-search-time.txt"
    pro = "MEB"
    dataset = ["writer", "WikiPedia"]
    title = ["Writers", "WikiPedia"]

    for i in range(0, len(dataset)):
        da = dataset[i]
        arrayY = []
        maxy = 0
        for j in range(0, len(algo)):
            al = algo[j]
            if j <= 1:
                path = path_baseline
            else:
                path = path_algo
            X, Y, Csize_array, up = read_data2(path, da, al, pro)
            arrayY.append(Y)
            maxy = max(maxy, up)
        print("dataset:", da)
        array = np.transpose(arrayY)  # array.T
        Draw5(array, tu_lab, title[i])

    # right 0.76
    # top 0.83
    # cut 1949x1195

