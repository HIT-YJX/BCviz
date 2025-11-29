import scipy.stats as stats

path = ""
NUM = 500000  # for id in row 1


def BCviz_similarity(files):
    X_array = []
    for dataset in files:
        dataset = path + dataset
        file = open(dataset, "r")
        id = []
        for i, d in enumerate(file):
            datalist = d.strip().split()
            if len(datalist) < 2:
                continue
            if int(datalist[0]) == 0:
                id.append(int(datalist[1]))
            else:
                id.append(int(datalist[1])+NUM)
            X_array.append(id)

    # print(kendalltau(X_array[0], X_array[2]))
    result = stats.kendalltau(X_array[0], X_array[2])
    print("the similarity between BCviz+ and BCviz: ", result.correlation)

    # print(kendalltau(X_array[1], X_array[2]))
    result = stats.kendalltau(X_array[1], X_array[2])
    print("the similarity between BCviz- and BCviz: ", result.correlation)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print("---Index Similarity in terms of Kendall's tau----")
    file_tail = ["_MEB_BCviz+.txt", "_MEB_BCviz-.txt", "_MEB_BCviz.txt"]
    data_name = input("Input the dataset name:")  # e.g. marvel, WikiPedia
    path = input("Input the absolute path of the index directory:")  # e.g. D:/Source Code/BCviz/Index-results/

    file_name = []
    for i in range(0, len(file_tail)):
        file_name.append(data_name + file_tail[i])
    print("dataset:", data_name)
    BCviz_similarity(file_name)
