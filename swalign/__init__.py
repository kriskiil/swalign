#!/usr/bin/env python
'''
Simple Smith-Waterman aligner
'''
import sys, StringIO, bisect, math

class ScoringMatrix(object):
    '''
    Read scoring matrix from a file or string

    Matrix should be space-delimited in a format like:

      A C G T
    A 1 0 0 0
    C 0 1 0 0
    G 0 0 1 0
    T 0 0 0 1

    Rows and Columns must be in the same order

    '''
    def __init__(self, filename=None, text=None, wildcard_score=0):
        assert filename or text

        if filename:
            fs = open(filename)
        else:
            fs = StringIO.StringIO(text)

        self.scores = []
        self.bases = None
        self.wildcard_score = wildcard_score

        for line in fs:
            if line[0] == '#' or not line.strip():
                continue

            if not self.bases:
                self.bases = line.split()
                self.base_count = len(self.bases)
            else:
                cols = line.split()
                self.scores.extend([float(x) for x in cols[1:]])

        fs.close()

    def score(self, one, two, wildcard=None):
        if self.wildcard_score and wildcard and (one in wildcard or two in wildcard):
            return self.wildcard_score

        one_idx = 0
        two_idx = 0
        for i, b in enumerate(self.bases):
            if b == one:
                one_idx = i
            if b == two:
                two_idx = i

        return self.scores[(one_idx * self.base_count) + two_idx]

class NumPyMatrix(object):
    def __init__(self, rows, cols, init=None):
        self.rows = rows
        self.cols = cols

        self.a = array(full((rows, cols), init[0], dtype=np.int))
        self.b = array(full((rows, cols), init[1], dtype=np.str))
        self.c = array(full((rows, cols), init[2], dtype=np.int))

    def get(self, row, col):
        return (self.a[row,col], self.b[row,col], self.c[row,col])

    def set(self, row, col, val):
        self.a[row,col] = val[0]
        self.b[row,col] = val[1]
        self.c[row,col] = val[2]

class IdentityScoringMatrix(object):
    def __init__(self, match=1, mismatch=-1):
        self.match = match
        self.mismatch = mismatch

    def score(self, one, two, wildcard=None):
        if wildcard and (one in wildcard or two in wildcard):
            return self.match

        if one == two:
            return self.match
        return self.mismatch

NucleotideScoringMatrix = IdentityScoringMatrix

class Matrix(object):
    def __init__(self, rows, cols, init=None):
        self.rows = rows
        self.cols = cols
        self.values = [init, ] * rows * cols

    def get(self, row, col):
        return self.values[(row * self.cols) + col]

    def set(self, row, col, val):
        self.values[(row * self.cols) + col] = val

class AlignmentPool:
    """
    Kind of like a priority queue of tuples except once it fills up, the element
    with the lowest "priority" will be discarded
    """
    def __init__(self, max_size, threshold):
        self.max_size = max_size
        self.threshold = threshold
        self.valid_count = 0
        self.list = []
    def put(self, elem):
        if elem[0] >= self.threshold:
            self.valid_count += 1
            #ordering cheat
            elem = self._change_ordering(elem)
            #insert in order
            bisect.insort_left(self.list, elem)
            #if we overflow
            if(len(self.list) > self.max_size):
                del self.list[self.max_size-1] #delete the overflowing elem
    def get(self):
        try:
            elem = self._change_ordering(self.list[0])
            #remove top element
            self.list = self.list[1:]
            return elem
        except IndexError, e: #not workign!
            return None
        #restore ordering values
    def empty(self):
        return not self.list
    def _change_ordering(self, elem):
        elem = list(elem)
        elem[0] *= -1
        return tuple(elem)

class LocalAlignment(object):
    def __init__(self, scoring_matrix, gap_penalty=-1, num_alignments=10, alignment_quality=90, gap_extension_penalty=-1, gap_extension_decay=0.0, prefer_gap_runs=True, verbose=False, globalalign=False, wildcard=None, full_query=False):
        self.scoring_matrix = scoring_matrix
        self.gap_penalty = gap_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.gap_extension_decay = gap_extension_decay
        self.verbose = verbose
        self.prefer_gap_runs = prefer_gap_runs
        self.globalalign = globalalign
        self.wildcard = wildcard
        self.full_query = full_query
        self.num_alignments = num_alignments
        self.alignment_quality = alignment_quality

    def align(self, ref, query, ref_name='', query_name='', rc=False):
        #calculate acceptance threshold for countable alignment
        threshold = math.floor(((len(query)*self.scoring_matrix.match)/100.0)*self.alignment_quality)
        #pool to hold the top num_alignments alignments
        pool = AlignmentPool(self.num_alignments, threshold)
        orig_ref = ref
        orig_query = query

        ref = ref.upper()
        query = query.upper()

        matrix = Matrix(len(query) + 1, len(ref) + 1, (0, ' ', 0))
        for row in xrange(1, matrix.rows):
            matrix.set(row, 0, (0, 'i', 0))

        for col in xrange(1, matrix.cols):
            matrix.set(0, col, (0, 'd', 0))

        max_val = 0
        max_row = 0
        max_col = 0

        # calculate matrix
        for row in xrange(1, matrix.rows):
            for col in xrange(1, matrix.cols):
                mm_val = matrix.get(row - 1, col - 1)[0] + self.scoring_matrix.score(query[row - 1], ref[col - 1], self.wildcard)

                ins_run = 0
                del_run = 0

                if matrix.get(row - 1, col)[1] == 'i':
                    ins_run = matrix.get(row - 1, col)[2]
                    if matrix.get(row - 1, col)[0] == 0:
                        # no penalty to start the alignment
                        ins_val = 0
                    else:
                        if not self.gap_extension_decay:
                            ins_val = matrix.get(row - 1, col)[0] + self.gap_extension_penalty
                        else:
                            ins_val = matrix.get(row - 1, col)[0] + min(0, self.gap_extension_penalty + ins_run * self.gap_extension_decay)
                else:
                    ins_val = matrix.get(row - 1, col)[0] + self.gap_penalty

                if matrix.get(row, col - 1)[1] == 'd':
                    del_run = matrix.get(row, col - 1)[2]
                    if matrix.get(row, col - 1)[0] == 0:
                        # no penalty to start the alignment
                        del_val = 0
                    else:
                        if not self.gap_extension_decay:
                            del_val = matrix.get(row, col - 1)[0] + self.gap_extension_penalty
                        else:
                            del_val = matrix.get(row, col - 1)[0] + min(0, self.gap_extension_penalty + del_run * self.gap_extension_decay)

                else:
                    del_val = matrix.get(row, col - 1)[0] + self.gap_penalty

                if self.globalalign or self.full_query:
                    cell_val = max(mm_val, del_val, ins_val)
                else:
                    cell_val = max(mm_val, del_val, ins_val, 0)

                if not self.prefer_gap_runs:
                    ins_run = 0
                    del_run = 0

                if del_run and cell_val == del_val:
                    val = (cell_val, 'd', del_run + 1)
                elif ins_run and cell_val == ins_val:
                    val = (cell_val, 'i', ins_run + 1)
                elif cell_val == mm_val:
                    val = (cell_val, 'm', 0)
                elif cell_val == del_val:
                    val = (cell_val, 'd', 1)
                elif cell_val == ins_val:
                    val = (cell_val, 'i', 1)
                else:
                    val = (0, 'x', 0)

                if val[0] >= max_val and val[0] > 0:
                    max_val = val[0]
                    max_row = row
                    max_col = col
                    #add to the pool
                    pool.put((max_val, max_row, max_col))

                matrix.set(row, col, val)
        #list for final alignments
        alignments = []
        #process each alignment in the pool
        while not pool.empty():
            current = pool.get()
            max_val = current[0]
            max_row = current[1]
            max_col = current[2]
            # backtrack
            if self.globalalign:
                # backtrack from last cell
                row = matrix.rows - 1
                col = matrix.cols - 1
                val = matrix.get(row, col)[0]
            elif self.full_query:
                # backtrack from max in last row
                row = matrix.rows - 1
                max_val = 0
                col = 0
                for c in xrange(1, matrix.cols):
                    if matrix.get(row, c)[0] > max_val:
                        col = c
                        max_val = matrix.get(row, c)[0]
                col = matrix.cols - 1
                val = matrix.get(row, col)[0]
            else:
                # backtrack from max
                row = max_row
                col = max_col
                val = max_val

            op = ''
            aln = []

            path = []
            while True:
                val, op, runlen = matrix.get(row, col)

                if self.globalalign:
                    if row == 0 and col == 0:
                        break
                elif self.full_query:
                    if row == 0:
                        break
                else:
                    if val <= 0:
                        break

                path.append((row, col))
                aln.append(op)

                if op == 'm':
                    row -= 1
                    col -= 1
                elif op == 'i':
                    row -= 1
                elif op == 'd':
                    col -= 1
                else:
                    break

            aln.reverse()
            if self.verbose:
                self.dump_matrix(ref, query, matrix, path)
                print aln
                print (max_row, max_col), max_val

            cigar = _reduce_cigar(aln)

            alignment = Alignment(orig_query, orig_ref, row, col, cigar, max_val, ref_name, query_name, rc, self.globalalign, self.wildcard)
            alignments.append(alignment)

        return (pool.valid_count, alignments)

    def dump_matrix(self, ref, query, matrix, path, show_row=-1, show_col=-1):
        print('      -      ')
        print('       '.join(ref))
        print('\n')
        for row in xrange(matrix.rows):
            if row == 0:
                print('-')
            else:
                print(query[row - 1])

            for col in xrange(matrix.cols):
                if show_row == row and show_col == col:
                    print('       *')
                else:
                    print(' %5s%s%s' % (matrix.get(row, col)[0], matrix.get(row, col)[1], '$' if (row, col) in path else ' '))
            print('\n')

def _reduce_cigar(operations):
    count = 1
    last = None
    ret = []
    for op in operations:
        if last and op == last:
            count += 1
        elif last:
            ret.append((count, last.upper()))
            count = 1
        last = op

    if last:
        ret.append((count, last.upper()))
    return ret

def _cigar_str(cigar):
    out = ''
    for num, op in cigar:
        out += '%s%s' % (num, op)
    return out

class Alignment(object):
    def __init__(self, query, ref, q_pos, r_pos, cigar, score, ref_name='', query_name='', rc=False, globalalign=False, wildcard=None):
        self.compressed = False
        self.q_pos = q_pos
        self.r_pos = r_pos
        self.cigar = cigar
        self.score = score
        self.r_name = ref_name
        self.q_name = query_name
        self.rc = rc
        self.globalalign = globalalign
        self.wildcard = wildcard

        self.r_offset = 0
        self.r_region = None

        self.orig_query = query
        self.query = query.upper()

        self.orig_ref = ref
        self.ref = ref.upper()

        q_len = 0
        r_len = 0

        self.matches = 0
        self.mismatches = 0

        i = self.r_pos
        j = self.q_pos

        for count, op in self.cigar:
            if op == 'M':
                q_len += count
                r_len += count
                for k in xrange(count):
                    if self.query[j] == self.ref[i]:
                        self.matches += 1
                    else:
                        self.mismatches += 1
                    i += 1
                    j += 1

            elif op == 'I':
                q_len += count
                j += count
                self.mismatches += count
            elif op == 'D':
                r_len += count
                i += count
                self.mismatches += count

        self.q_end = q_pos + q_len
        self.r_end = r_pos + r_len
        if self.mismatches + self.matches > 0:
            self.identity = float(self.matches) / (self.mismatches + self.matches)
        else:
            self.identity = 0

    # def compress(self):
    #     """
    #     Deletes the objects copy of the genomes to save memory to allow lots
    #     of alignments to be stored in a list in memory
    #     """
    #     self.ref = None
    #     self.orig_ref = None
    #     self.query = None
    #     self.orig_quer = None
    #     self.compressed = True
    #
    # def decompress(self, query, ref):
    #     """
    #     Adds back the genomes being aligned
    #     """
    #     self.orig_query = query
    #     self.query = query.upper()
    #
    #     self.orig_ref = ref
    #     self.ref = ref.upper()
    #     self.compressed = False

    def set_ref_offset(self, ref, offset, region):
        self.r_name = ref
        self.r_offset = offset
        self.r_region = region

    @property
    def extended_cigar_str(self):
        if self.compressed:
            raise Exception('Alignment is compressed')
        qpos = 0
        rpos = 0
        ext_cigar_str = ''
        working = []
        for count, op in self.cigar:
            if op == 'M':
                for k in xrange(count):
                    if self.query[self.q_pos + qpos + k] == self.ref[self.r_pos + rpos + k]:
                        ext_cigar_str += 'M'
                    else:
                        ext_cigar_str += 'X'
                qpos += count
                rpos += count

            elif op == 'I':
                qpos += count
                ext_cigar_str += 'I' * count
            elif op == 'D':
                rpos += count
                ext_cigar_str += 'D' * count

            working = _reduce_cigar(ext_cigar_str)

        out = ''
        for num, op in working:
            out += '%s%s' % (num, op)
        return out

    @property
    def cigar_str(self):
        if self.compressed:
            raise Exception('Alignment is compressed')
        return _cigar_str(self.cigar)

    def dump(self, wrap=None, out=sys.stdout):
        if self.compressed:
            raise Exception('Alignment is compressed')
        i = self.r_pos
        j = self.q_pos

        q = ''
        m = ''
        r = ''
        qlen = 0
        rlen = 0

        for count, op in self.cigar:
            if op == 'M':
                qlen += count
                rlen += count
                for k in xrange(count):
                    q += self.orig_query[j]
                    r += self.orig_ref[i]
                    if self.query[j] == self.ref[i] or (self.wildcard and (self.query[j] in self.wildcard or self.ref[i] in self.wildcard)):
                        m += '|'
                    else:
                        m += '.'

                    i += 1
                    j += 1
            elif op == 'D':
                rlen += count
                for k in xrange(count):
                    q += '-'
                    r += self.orig_ref[i]
                    m += ' '
                    i += 1
            elif op == 'I':
                qlen += count
                for k in xrange(count):
                    q += self.orig_query[j]
                    r += '-'
                    m += ' '
                    j += 1

            elif op == 'N':
                q += '-//-'
                r += '-//-'
                m += '    '

        if self.q_name:
            out.write('Query: %s%s (%s nt)\n' % (self.q_name, ' (reverse-compliment)' if self.rc else '', len(self.query)))
        if self.r_name:
            if self.r_region:
                out.write('Ref  : %s (%s)\n\n' % (self.r_name, self.r_region))
            else:
                out.write('Ref  : %s (%s nt)\n\n' % (self.r_name, len(self.ref)))

        poslens = [self.q_pos + 1, self.q_end + 1, self.r_pos + self.r_offset + 1, self.r_end + self.r_offset + 1]
        maxlen = max([len(str(x)) for x in poslens])
        output = {}
        q_pre = 'Query: %%%ss ' % maxlen
        r_pre = 'Ref  : %%%ss ' % maxlen
        m_pre = ' ' * (8 + maxlen)
        rpos = self.r_pos
        if not self.rc:
            qpos = self.q_pos
        else:
            qpos = self.q_end
        output['qs'] = qpos
        while q and r and m:
            if not self.rc:
                out.write(q_pre % (qpos + 1))  # pos is displayed as 1-based
            else:
                out.write(q_pre % (qpos))  # revcomp is 1-based on the 3' end

            if wrap:
                qfragment = q[:wrap]
                mfragment = m[:wrap]
                rfragment = r[:wrap]

                q = q[wrap:]
                m = m[wrap:]
                r = r[wrap:]
            else:
                qfragment = q
                mfragment = m
                rfragment = r

                q = ''
                m = ''
                r = ''

            out.write(qfragment)
            if not self.rc:
                for base in qfragment:
                    if base != '-':
                        qpos += 1
            else:
                for base in qfragment:
                    if base != '-':
                        qpos -= 1

            if not self.rc:
                out.write(' %s\n' % qpos)
            else:
                out.write(' %s\n' % (qpos + 1))
            output['qe'] = qpos
            out.write(m_pre)
            out.write(mfragment)
            out.write('\n')
            out.write(r_pre % (rpos + self.r_offset + 1))
            out.write(rfragment)
            output['rs'] = rpos
            for base in rfragment:
                if base != '-':
                    rpos += 1
            output['re'] = rpos
            output['c'] = self.cigar_str
            out.write(' %s\n\n' % (rpos + self.r_offset))
        output['s'] = self.score
        out.write("Score: %s\n" % self.score)
        out.write("Matches: %s (%.1f%%)\n" % (self.matches, self.identity * 100))
        output['m'] = self.matches
        out.write("Mismatches: %s\n" % (self.mismatches,))
        output['mm'] = self.mismatches
        out.write("CIGAR: %s\n" % self.cigar_str)
        return output
def fasta_gen(fname):
    def gen():
        seq = ''
        name = ''
        comments = ''

        if fname == '-':
            f = sys.stdin
            name = 'stdin'
        else:
            f = open(fname)

        for line in f:
            if line[0] == '>':
                if name and seq:
                    yield (name, seq, comments)

                spl = line[1:].strip().split(' ', 1)
                name = spl[0]
                if len(spl) > 1:
                    comments = spl[1]
                else:
                    comments = ''

                seq = ''
            else:
                seq += line.strip()

        if name and seq:
            yield (name, seq, comments)

        if fname != '-':
            f.close()
    return gen

def seq_gen(name, seq):
    def gen():
        yield (name, seq, '')

    return gen

def extract_region(comments):
    ref = None
    start = None
    # start_offset = 0
    # end_offset = 0

    try:
        attrs = comments.split(' ')
        for attr in attrs:
            if '=' in attr:
                k, v = attr.split('=')
                if k == 'range':
                    spl = v.split(':')
                    ref = spl[0]
                    start, end = [int(x) for x in spl[1].split('-')]
                # elif k == "5'pad":
                #     start_offset = int(v)
                # elif k == "3'pad":
                #     end_offset = int(v)
    except:
        pass

    if ref and start:
        return (ref, start - 1, '%s:%s-%s' % (ref, start, end))

    return None

__revcomp = {}
for a, b in zip('atcgATCGNn', 'tagcTAGCNn'):
    __revcomp[a] = b
__cache = {}


def revcomp(seq):
    if seq in __cache:
        return __cache[seq]

    ret = []
    for s in seq.upper()[::-1]:
        ret.append(__revcomp[s])

    __cache[seq] = ''.join(ret)
    return __cache[seq]


#     sw.align('ACACACTA','AGCACACA').dump()
#     aln=sw.align("AAGGGGAGGACGATGCGGATGTTC","AGGGAGGACGATGCGG")
#     aln.dump()
