#Tryptic peptide generator
#Logic credit: Santosh Kumar Behera, beheras40@gmail.com
#Code credit: Sandeep Kasaragod, sandeep.kolya@gmail.com

class ProteinDigestion(object):

    def __init__(self,sequence,missed_clevage,min_len,max_len):
        self.sequence = sequence
        self.missed_clevage = missed_clevage
        self.min_len = min_len
        self.max_len = max_len
        
    def trypsin(self):
        if 'K' in self.sequence or 'R' in self.sequence:
            get_dup_k = [i for i in range(len(self.sequence)) if self.sequence.startswith('K', i)]
            get_dup_r = [j for j in range(len(self.sequence)) if self.sequence.startswith('R', j)]
            merge_list = sorted(get_dup_k + get_dup_r)
            merge_list_fltrd = [i for i in merge_list if i+1 < len(self.sequence) and self.sequence[i + 1] !='P'] #look for KP or RP position
            merge_list_fltrd.append(len(self.sequence))
            initialize = 0
            for iter_lst in range(len(merge_list_fltrd) - int(self.missed_clevage)):
                peptide = (self.sequence[initialize: int(merge_list_fltrd[iter_lst + self.missed_clevage]) + 1])
                if len(peptide) >= int(self.min_len) and len(peptide) <= int(self.max_len):
                    yield peptide
                initialize = merge_list_fltrd[iter_lst] + 1

    def lysc(self):
        get_dup_k = [i for i in range(len(self.sequence)) if self.sequence.startswith('K', i)]
        merge_list_fltrd = get_dup_k
        merge_list_fltrd.append(len(self.sequence))
        initialize = 0
        for iter_lst in range(len(merge_list_fltrd) - int(self.missed_clevage)):
            peptide = (self.sequence[initialize: int(merge_list_fltrd[iter_lst + self.missed_clevage]) + 1])
            if len(peptide) >= int(self.min_len) and len(self.peptide) <= int(self.max_len):
                yield peptide
            initialize = merge_list_fltrd[iter_lst] + 1
            
    def chymotrypsin(self):
        get_dup_f = [i for i in range(len(self.sequence)) if self.sequence.startswith('F', i)]
        get_dup_w = [j for j in range(len(self.sequence)) if self.sequence.startswith('W', j)]
        get_dup_y = [j for j in range(len(self.sequence)) if self.sequence.startswith('Y', j)]
        merge_list = sorted(get_dup_f + get_dup_w + get_dup_y)
        merge_list_fltrd = [i for i in merge_list if i+1 < len(self.sequence) and self.sequence[i + 1] !='P'] #look for KP or RP position
        merge_list_fltrd.append(len(self.sequence))
        initialize = 0
        for iter_lst in range(len(merge_list_fltrd) - int(self.missed_clevage)):
            peptide = (self.sequence[initialize: int(merge_list_fltrd[iter_lst + self.missed_clevage]) + 1])
            if len(peptide) >= int(self.min_len) and len(self.peptide) <= int(self.max_len):
                yield peptide
            initialize = merge_list_fltrd[iter_lst] + 1


class decoy_pep(object):
    def __init__(self, peptide):
        self.peptide = peptide

    def Reverse(self):
        pep = [self.peptide[p] for p in range(len(self.peptide)) if p+1 == len(self.peptide)] + [self.peptide[p] for p in range(len(self.peptide)) if p+1 != len(self.peptide)]
        pep.reverse()
        return ''.join(pep)
        
