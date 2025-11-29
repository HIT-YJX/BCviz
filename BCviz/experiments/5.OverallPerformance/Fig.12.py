import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)
sns.set_style("dark")
sns.set_style("ticks")
colorlist = ["#006662", "#006694", "#7a5195", "#ef5675", "#ffa600"]  # light

###  index  ###
# construct_BCviz = []
algo_index = ["BCviz", "BCviz+", "BCviz-"]
algo_lab_index = ["$MBCI$", "$BCviz$", "$BCviz+$", "$BCviz-$"]

###  search  ###
algo_baseline = [["MBC", "MBCI"], ["BIIP", "MBC"], ["ooMBC", "MBC"]]
algo = ["BCviz", "BCviz+", "BCviz-"]
algo_lab = [["$MBC$", "$MBCI$", "$BCviz$", "$BCviz+$", "$BCviz-$"],
            ["$BIIP$", "$MBC$", "$BCviz$", "$BCviz+$", "$BCviz-$"],
            ["$ooMBB$", "$MBC$", "$BCviz$", "$BCviz+$", "$BCviz-$"]]

datasets = ["marvel", "writer", "BookCrossing",
            "Team", "Actor-movie", "Twitter",
            "WikiPedia", "DBLP"]
datasets_lab = ["$MV$", "$WT$", "$BC$",
            "$TM$", "$AM$", "$TT$",
            "$WP$", "$DB$"]


def read_data_construct(path, sizetype, algo_name):
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
                maxy = max(maxy, float(datalist[5]))  # construct
                y.append(float(datalist[5]))
                x.append(row)
                # Csize.append(int(datalist[5]))
                row = row + 1
    return x, y, Csize, maxy


def read_data(path, pro_type, algo_name):
    file = open(path, "r")
    x = []
    y = []
    Csize = []
    row = 0
    maxy = 0
    for i, d in enumerate(file):
        datalist = d.strip().split()
        if len(datalist) < 4:
            continue
        if datalist[0] == pro_type and datalist[1] == algo_name:
            if float(datalist[3]) >= 0.1:
                maxy = max(maxy, float(datalist[3]))
                y.append(float(datalist[3]))
            else:
                y.append(0.1)
            x.append(row)
            # Csize.append(int(datalist[3]))
            row = row + 1
    return x, y, Csize, maxy


def Draw5(rcount, x_label, title, Buttom, algo_labels):
    plt.figure(figsize=(25, 13))
    plt.rcParams['figure.figsize'] = (80, 60)
    plt.xlabel("", fontsize=60)
    plt.xticks(fontsize=48)  # 28 40
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=True
    )
    plt.ylabel("Total Time(s)", fontsize=60)  # 49 80
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

    if pro_type == "MEB":
        s2 = [i[0] for i in Buttom]
        s3 = [i[1] for i in Buttom]
        s4 = [i[2] for i in Buttom]
        s5 = [i[3] for i in Buttom]
    else:
        s3 = [i[0] for i in Buttom]
        s4 = [i[1] for i in Buttom]
        s5 = [i[2] for i in Buttom]

    y1 = [i[0] for i in rcount]
    y2 = [i[1] for i in rcount]
    y3 = [i[2] for i in rcount]
    y4 = [i[3] for i in rcount]
    y0 = [0, 0, 0, 0, 0, 0, 0, 0]
    y5 = [i[4] for i in rcount]

    b0 = plt.bar(x1, y0, width=rwidth, color='white', alpha=0.0, label='Search', log=1)
    b1 = plt.bar(x1, y1, width=rwidth, color=colorlist[1], lw=0.5, edgecolor='black', alpha=0.9, label=algo_labels[0], log=1)
    b2 = plt.bar(x2, y2, width=rwidth, color=colorlist[0], lw=0.5, edgecolor='black', alpha=0.9, label=algo_labels[1], log=1)
    b3 = plt.bar(x3, y3, width=rwidth, color=colorlist[2], lw=0.5, edgecolor='black', alpha=0.9, label=algo_labels[2], log=1, tick_label=x_label)
    b4 = plt.bar(x4, y4, width=rwidth, color=colorlist[3], lw=0.5, edgecolor='black', alpha=0.9, label=algo_labels[3], log=1)
    b5 = plt.bar(x5, y5, width=rwidth, color=colorlist[4], lw=0.5, edgecolor='black', alpha=0.9, label=algo_labels[4], log=1)

    t0 = plt.bar(x2, y0, width=rwidth, color='white', alpha=0.0, label='Construction', log=1)
    if title == "MEB":
        t2 = plt.bar(x2, s2, bottom=y2, width=rwidth, color=colorlist[0], hatch='\\\\\\', lw=0.8, edgecolor='black', alpha=0.3, label=algo_labels[1], log=1)
    else:
        t2 = plt.bar(x2, y0, width=rwidth, color='white', alpha=0.0, label='.', log=1)
    t3 = plt.bar(x3, s3, bottom=y3, width=rwidth, color=colorlist[2], hatch='\\\\\\', lw=0.8, edgecolor='black', alpha=0.3, label=algo_labels[2], log=1, tick_label=x_label)
    t4 = plt.bar(x4, s4, bottom=y4, width=rwidth, color=colorlist[3], hatch='\\\\\\', lw=0.8, edgecolor='black', alpha=0.3, label=algo_labels[3], log=1)
    t5 = plt.bar(x5, s5, bottom=y5, width=rwidth, color=colorlist[4], hatch='\\\\\\', lw=0.8, edgecolor='black', alpha=0.3, label=algo_labels[4], log=1)

    ax = plt.gca()
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    bwith = 2  # 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    ax.set_yticks([0.1, 1, 10, 100, 1000, 10000, 100000])
    plt.yticks([0.1, 1, 10, 100, 1000, 10000, 100000], ['$10^{-1}$', '$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$', '$10^{4}$', '$10^{5}$'])

    if title == "MEB":
        plt.legend(handles=[b0, b1, b2, b3, b4, b5, t0, t2, t3, t4, t5], loc='best', frameon=True, fontsize=45, ncol=2) # , bbox_to_anchor=(0.45, 0.5)
    else:
        plt.legend(handles=[b0, b1, b2, b3, b4, b5, t0, t3, t4, t5, t2], loc='best', frameon=True, fontsize=45, ncol=2) # , bbox_to_anchor=(0.45, 0.5)

    plt.xticks()
    file = title + "-total-time.jpg"
    plt.savefig(file)
    print(file + " has been generated.")
    # plt.show()


if __name__ == '__main__':
    print("--------Overall Performance--------")
    path_index = "index-time.txt"
    path_search = "total-time.txt"
    path_search_baseline = "baseline-total-time.txt"
    problems = ["MEB", "MVB", "MBB"]

    for i in range(0, len(problems)):
        pro_type = problems[i]
        print("Problem type:", pro_type)

        # once index time
        arrayButtom = []
        if pro_type == "MEB":
            X, Y, Csize_array, up = read_data_construct("baseline-index-time.txt", pro_type, "MBCI")
            arrayButtom.append(Y)
        for al in algo_index:
            X, Y, Csize_array, up = read_data_construct(path_index, pro_type, al)
            arrayButtom.append(Y)
        arrayB = np.transpose(arrayButtom)  # array.T
        # print("index time:")
        # print(arrayB)

        # multiple search time
        arrayY = []
        for al in algo_baseline[i]:
            X, Y, Csize_array, up = read_data(path_search_baseline, pro_type, al)
            arrayY.append(Y)
        for al in algo:
            X, Y, Csize_array, up = read_data(path_search, pro_type, al)
            arrayY.append(Y)
        array = np.transpose(arrayY)  # array.T
        # print("total search time:")
        # print(array)
        Draw5(array, datasets_lab, pro_type, arrayB, algo_lab[i])

    # tall 0.7 right 0.998
    # 2308x1010
