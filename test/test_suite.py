from proteins import make_amino_dicts, make_one_hot_res

def test_one_hot1():
    sequence='ACDEF'
    class FakeInvariant:
        def __init__(self, sequence='ACDEF'):
            self.sequence = sequence
            self.length = len(sequence)

    res_to_idx, _  = make_amino_dicts()

    one_protein = FakeInvariant(sequence)
    one_hot = make_one_hot_res(one_protein, res_to_idx)
    assert one_hot.shape[0] == len(sequence)
    assert one_hot.shape[1] == len(res_to_idx)
    for i in range(len(sequence)):
        assert one_hot[i,i] == 1
    for i in range(one_hot.shape[0]):
        for j in range(one_hot.shape[1]):
            if i<len(sequence) and i==j:
                assert one_hot[i,j] ==1
            else:
                assert one_hot[i,j]==0

def test_one_hot2():
    sequence='YAAAAAAAAAAA'
    class FakeInvariant:
        def __init__(self, sequence=''):
            self.sequence = sequence
            self.length = len(sequence)

    res_to_idx, _  = make_amino_dicts()

    one_protein = FakeInvariant(sequence)
    one_hot = make_one_hot_res(one_protein, res_to_idx)
    assert one_hot.shape[0] == len(sequence)
    assert one_hot.shape[1] == len(res_to_idx)

    for i in range(one_hot.shape[0]):
        for j in range(one_hot.shape[1]):
            if i==0:
                if j==len(res_to_idx)-1:
                    assert one_hot[i,j] ==1
                else:
                    assert one_hot[i,j] ==0
            elif j==0:
                assert one_hot[i,j]==1
            else:
                assert one_hot[i,j]==0