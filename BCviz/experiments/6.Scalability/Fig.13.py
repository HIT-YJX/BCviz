import matplotlib.pyplot as plt
import math


maker = ['d', '*', 'o', 's']
colorlist = ["#006662", "#7a5195", "#ef5675", "#ffa600"]  # light
rate = ["20", "40", "60", "80", ""]
rate_label = ["20%", "40%", "60%", "80%", "100%"]


# def exponential_func(x):
#     y = math.log10(x)
#     return y


def read_data(path, da, p1, algo_name, p2, pv):
    file = open(path, "r")
    x = []
    y = []
    Csize = []
    row = 0
    maxy = 0

    r = 0
    data_name = da + rate[r]
    for i, d in enumerate(file):
        datalist = d.strip().split()
        if len(datalist) < pv+1:
            continue
        if datalist[p1] == data_name:
            if datalist[p2] == algo_name:
                value = float(datalist[pv])
                maxy = max(maxy, value)  # construct
                y.append(value)
                x.append(row)
                r = r + 1
                if r > 4:
                    break
                data_name = da + rate[r]
        # Csize.append(int(datalist[5]))
        row = row + 1
    return x, y, Csize, maxy


def draw_line(X, Y, labs, title):
    plt.figure(figsize=(30, 19))
    plt.rcParams['figure.figsize'] = (50, 30)
    rwidth = 20
    x = [i*rwidth+rwidth for i in X]
    y0 = ["" for i in X]

    plt.title(title, size=70)
    plt.xticks(fontsize=70)
    plt.yticks(fontsize=70)
    plt.tick_params(axis='y', width=4, length=30)

    maxy = 0
    miny = 100000
    ax2 = plt.subplot(1, 1, 1)  # , facecolor="ghostwhite"
    for i in range(0, len(Y)):
        ax2.plot(x, Y[i], lw=6, ls='-', c=colorlist[i], marker=maker[i], markersize=30, label=labs[i])
        for j in range(0, len(Y[i])):
            maxy = max(maxy, Y[i][j])
            miny = min(miny, Y[i][j])
    ax2.set_xticks(x)
    plt.xticks(x, rate_label)

    ax2.legend(loc='best', fontsize=60)

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

    ax.set_yscale('log')
    ax2.set_xlabel('Edge Percentage', fontsize=70)  # 18
    ax2.set_ylabel('Construction Time(s)', fontsize=70)  # 18
    file = title + "-scalability.jpg"
    plt.savefig(file)
    print(file + " has been generated.")
    # plt.show()


if __name__ == '__main__':
    print("---Scalability (Index Construction Time)---")
    path = "scalability.txt"
    path_baseline = "baseline-scalability.txt"
    dataset = ["marvel", "Twitter"]
    dataset_lab = ["Marvel", "Twitter"]
    algo = ["MBCI", "BCviz", "BCviz+", "BCviz-"]
    algo_lab = ["$MBCI$", "$BCviz$", "$BCviz+$", "$BCviz-$"]

    row = 5
    for i in range(0, len(dataset)):
        da = dataset[i]
        Y = []
        print("dataset:", da)

        x, y, Csize_array, up = read_data(path_baseline, da, 2, algo[0], 1, row)
        Y.append(y)

        x, y, Csize_array, up = read_data(path, da, 2, algo[1], 1, row)
        Y.append(y)
        x, y, Csize_array, up = read_data(path, da, 2, algo[2], 1, row)
        Y.append(y)
        x, y, Csize_array, up = read_data(path, da, 2, algo[3], 1, row)
        Y.append(y)

        draw_line(x, Y, algo_lab, dataset_lab[i])

    # save top 0.9 right 0.76
    # cut 2006x1302





