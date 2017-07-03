class Iterator:

    def __init__(self, depth, start, end):
        if len(start) < depth:
            start += [0]*(depth - len(start))
        if len(end) < depth:
            end += [0] * (depth - len(end))
        self.depth = depth
        self.start = list(start)
        self.end = list(end)
        self.index = list(start)
        self.start_flag = True

    def __next__(self):
        if self.start_flag:
            self.start_flag = False
            return self.index
        for i in range(self.depth):
            self.index[i] += 1
            if self.index[i] < self.end[i]:
                break
            else:
                self.index[i] = self.start[i]
                if i == self.depth - 1:
                    raise StopIteration
        return self.index

    def __iter__(self):
        return self


def main():
    list1 = [1,2]
    list2 = [4,6]
    depth = 2
    for item in Iterator(depth, list1, list2):
        print(item)

if __name__ == "__main__":
    main()