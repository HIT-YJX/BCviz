# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import matplotlib.pyplot as plt

colorlist = ["#ffa600", "#ef5675", "#7a5195"]  # light
NUM = 500000
path = ""


def BCviz(algo, files, name):
    # plt.rcParams['font.sans-serif'] = ['SimHei']  # show chinese
    plt.figure(figsize=(25, 13))
    plt.rcParams['figure.figsize'] = (80, 60)

    maxy = 0
    pos = 0
    X_array = []
    B = []
    for dataset in files:
        dataset = path + dataset
        file = open(dataset, "r")
        x = []
        id = []
        y = []
        row = 0
        my = 0
        for i, d in enumerate(file):
            datalist = d.strip().split()
            if len(datalist) < 2:
                continue
            my = max(my, int(datalist[2]))
            y.append(int(datalist[2]))
            if int(datalist[0]) == 0:
                id.append(int(datalist[1]))
            else:
                id.append(int(datalist[1])+NUM)
            X_array.append(id)
            x.append(row)
            row = row + 1
        # print(algo[pos], row)
        maxy = max(maxy, my)
        b = plt.plot(x, y, '-', color=colorlist[pos], alpha=0.8, linewidth=8, label=algo[pos])  # '#4169E1'
        B.append(b)
        plt.legend(loc="upper center", fontsize=55)
        pos = pos + 1

    ax = plt.gca()
    # set border color
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    # set border line-width
    bwith = 3
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    plt.ylim(0, maxy + maxy / 100)
    plt.xlabel('$vertex$', fontsize=70)
    plt.ylabel(r'$\sigma_\gamma(v)$', fontsize=70)
    plt.xticks(fontsize=55, visible=False)
    plt.yticks(fontsize=55)
    plt.title(name, size=70)

    file = name + "-IndexVisulization.jpg"
    plt.savefig(file)
    print(file + " has been generated.")
    # plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print("---------Index Visualization----------")
    datasets = ["marvel", "WikiPedia"]
    titles = ["Marvel", "WikiPedia"]
    file_tail = ["_MEB_BCviz+.txt", "_MEB_BCviz-.txt", "_MEB_BCviz.txt"]
    algo_lab = ["$BCviz+$", "$BCviz-$", "$BCviz$"]
    path = input("Input the absolute path of the index directory:")  # example, D:/Source Code/BCviz/Index-results/

    for j in range(0, len(datasets)):
        file_name = []
        label = []
        for i in range(0, len(file_tail)):
            file_name.append(datasets[j] + file_tail[i])
            label.append(algo_lab[i])

        print("dataset:", datasets[j])
        BCviz(label, file_name, titles[j])

    # top 0.8
    # 2376x1166 -> 2310x1156
