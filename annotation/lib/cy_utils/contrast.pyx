
import cython
cimport numpy as np
import numpy as np

@cython.profile(True)
@cython.cdivision(True)
cdef double calc_f1(double recall, double precision):
    """
    Static method to calculate the F1 statistic given precision
    and recall (order is unimportant). Definition:
    F1 = (2 * precision * recall) / (precision + recall)
    """
    cdef double result, summa, product

    if precision < 0 or recall < 0:
        raise ValueError("Negative values are an invalid input! ({0}, {1})".format(
            recall, precision))

    elif precision == 0 or recall == 0:
        return 0
    else:
        product = 2 * precision * recall
        summa = precision + recall
        result = product / summa
        return result


@cython.cdivision(True)
cpdef np.ndarray array_compare(np.ndarray[np.int_t] ref, np.ndarray[np.int_t] pred, double identity):

    # Calculate the junction intersection
    cdef:
        np.ndarray result
        char ccode
        np.ndarray ref_exons, pred_exons
        np.ndarray[np.int_t] ref_junc, pred_junc
        double junction_recall, junction_precision, junction_f1
        double exon_recall, exon_precision, exon_f1
        int cycling
        np.int_t junc_res
        np.ndarray exon_res
        double exon_found
        double splice_found

    ref_exons = ref.reshape((<int>ref.shape[0] / 2, 2))
    pred_exons = pred.reshape((<int>pred.shape[0] / 2, 2))

    found = 0
    for cycling in range(ref_exons.shape[0]):
        exon_res = np.where(np.all(ref_exons[cycling] == pred_exons, axis=1))[0]
        # print(res, res.shape[0])
        if exon_res.shape[0] > 0:
            found += 1

    exon_recall = found / < double > ref_exons.shape[0]
    exon_precision = found / < double > pred_exons.shape[0]
    exon_f1 = calc_f1(exon_recall, exon_precision)

    if ref.shape[0] > 2 and pred.shape[0] > 2:
        splice_found = 0
        ref_junc = ref[1:-1]
        pred_junc = pred[1:-1]
        splice_found = 0
        for cycling in range(ref_junc.shape[0]):
            junc_res = np.searchsorted(pred_junc, ref_junc[cycling])
            if junc_res >= pred_junc.shape[0]:
                break  # We have reached the end of comparisons
            if ref_junc[cycling] == pred_junc[junc_res]:
                splice_found += 1

        junction_recall = splice_found / <double> ref_junc.shape[0]
        junction_precision = splice_found / <double> pred_junc.shape[0]
        junction_f1 = calc_f1(junction_recall, junction_precision)
        if junction_f1 == 1:
            ccode = b'='
        elif junction_f1 > 0:
            ccode = b'j'
        else:
            ccode = b'h'
    else:
        if ref.shape[0] == 2 and pred.shape[0] == 2:
            junction_recall, junction_precision, junction_f1 = 1, 1, 1
            if identity >= 80:
                ccode = b'_'
            else:
                ccode = b'm'
        else:
            junction_recall, junction_precision, junction_f1 = 0, 0, 0
            if ref.shape[0] > pred.shape[0]:
                ccode = b'g'
            else:
                ccode = b'G'

    result = np.array([exon_recall, exon_precision, exon_f1,
                        junction_recall, junction_precision, junction_f1,
                       ccode])

    return result