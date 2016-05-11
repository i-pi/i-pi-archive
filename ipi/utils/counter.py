class Counter(object):

    def __init__(self):
        self.func_eval = 0
        
    def count(self):
        self.func_eval = self.func_eval + 1
        
counter = Counter()
counter1 = Counter()
