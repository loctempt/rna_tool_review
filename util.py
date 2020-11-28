class Util:
    @staticmethod
    def read_targets(target_file='targets.out'):
        target_dict = {}
        with open(target_file) as targets:
            lines = targets.readlines
            for line in lines:
                line = line.split()
                target_dict[line[0]] = line[1]
        return target_dict