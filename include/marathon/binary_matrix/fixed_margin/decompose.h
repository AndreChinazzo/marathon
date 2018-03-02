/*
 * Created on: May 12, 2017
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2017, Steffen Rechner
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#ifndef MARATHON_FIXED_MARGIN_DECOMPOSITION_H
#define MARATHON_FIXED_MARGIN_DECOMPOSITION_H

#include "marathon/binary_matrix/common.h"
#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/binary_matrix/fixed_margin/instance.h"
#include "marathon/binary_matrix/fixed_margin/realize.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Returns false, if sequence of row and column sums can be decomposed
             * into smaller sequences. Otherwise, return true.
             * @param rowsum Sequence of row sums.
             * @param colsum Sequence of column sums.
             * @param nrow Number of rows.
             * @param ncol Number of columns.
             * @return
             */
            inline
            bool isPrimitive(
                    const int *rowsum,
                    const int *colsum,
                    const int nrow,
                    const int ncol
            ) {

                // compute conjugated row sums
                int *rowsum_conjugated = new int[ncol];
                conjugate(rowsum_conjugated, rowsum, ncol, nrow);

                // descendingly sort column sums
                int *colsum_sorted = new int[ncol];
                memcpy(colsum_sorted, colsum, ncol * sizeof(int));
                std::sort(colsum_sorted, colsum_sorted + ncol, [](const int a, const int b) { return a > b; });

                bool reducible = false;

                // check dominance criterion
                int S1 = 0;
                int S2 = 0;
                for (int k = 0; k < ncol - 1; k++) {
                    S1 += rowsum_conjugated[k];
                    S2 += colsum_sorted[k];
                    if (S1 == S2) {
                        reducible = true;
                        break;
                    }
                }

                delete[] rowsum_conjugated;
                delete[] colsum_sorted;

                if (reducible)
                    return false;

                // compute conjugated column sums
                int *colsum_conjugated = new int[nrow];
                conjugate(colsum_conjugated, colsum, nrow, ncol);

                // descendingly sort row sums
                int *rowsum_sorted = new int[nrow];
                memcpy(rowsum_sorted, rowsum, nrow * sizeof(int));
                std::sort(rowsum_sorted, rowsum_sorted + nrow, [](const int a, const int b) { return a > b; });

                // check dominance criterion
                S1 = 0;
                S2 = 0;
                for (int k = 0; k < nrow - 1; k++) {
                    S1 += colsum_conjugated[k];
                    S2 += rowsum_sorted[k];
                    if (S1 == S2) {
                        reducible = true;
                        break;
                    }
                }

                delete[] colsum_conjugated;
                delete[] rowsum_sorted;

                return !reducible;
            }

            /**
             * Instance Decomposition Engine.
             */
            class Decompositor {

            protected:

                Instance _inst;
                bool _realizable;


                /**
                 * (Helper Function) Recursively decompose sequence seq.
                 * For each primitive sequence generated during the process,
                 * the function f is called.
                 * @param seq Realizable sequence pair.
                 * @param f Function called for each primitive sequence.
                 */
                void decompose_recursive(
                        const Instance &seq,
                        const std::function<void( const Instance &)> &f,
                        bool k1_may_exist,
                        bool k2_may_exist
                ) {

                    const int nrow = (int) seq.getNumRows();
                    const int ncol = (int) seq.getNumCols();

                    // std::cout << seq.toString() << std::endl;

                    // if seq may have a critical point k1
                    if (k1_may_exist) {

                        // try to find critical point k1

                        // build conjugate sequence
                        int *colsum_conjugate = new int[nrow];
                        conjugate(colsum_conjugate, &seq.colsum[0], nrow, ncol);

                        // determine critical point k1
                        int rowsum_sum = 0;
                        int colsum_conjugated_sum = 0;
                        int k1 = 0;
                        for (; k1 < nrow; k1++) {
                            rowsum_sum += seq.rowsum[k1];
                            colsum_conjugated_sum += colsum_conjugate[k1];
                            if (rowsum_sum == colsum_conjugated_sum) {
                                break;
                            }
                        }

                        // if critical point k1 does exist
                        if (k1 < nrow - 1) {

                            // sequence can be decomposed by splitting row indices at position k1
                            Instance sub1(k1 + 1, ncol);
                            Instance sub2(nrow - k1 - 1, ncol);
                            memcpy(&sub1.rowsum[0], &seq.rowsum[0], (k1 + 1) * sizeof(int));
                            memcpy(&sub2.rowsum[0], &seq.rowsum[k1 + 1], (nrow - k1 - 1) * sizeof(int));
                            memcpy(&sub1.rowindex[0], &seq.rowindex[0], (k1 + 1) * sizeof(int));
                            memcpy(&sub2.rowindex[0], &seq.rowindex[k1 + 1], (nrow - k1 - 1) * sizeof(int));
                            conjugate(&sub1.colsum[0], &colsum_conjugate[0], ncol, k1 + 1);
                            conjugate(&sub2.colsum[0], &colsum_conjugate[k1 + 1], ncol, nrow - k1 - 1);
                            memcpy(&sub1.colindex[0], &seq.colindex[0], ncol * sizeof(int));
                            memcpy(&sub2.colindex[0], &seq.colindex[0], ncol * sizeof(int));

                            decompose_recursive(sub1, f, false, true);
                            decompose_recursive(sub2, f, true, true);
                        } else { // critical point k1 does not exist

                            decompose_recursive(seq, f, false, k2_may_exist);
                        }

                        delete[] colsum_conjugate;

                    } else if (k2_may_exist) {

                        // try to find critical point k2

                        // build conjugate sequence
                        int *rowsum_conjugate = new int[ncol];
                        conjugate(rowsum_conjugate, &seq.rowsum[0], ncol, nrow);

                        // determine critical point k2
                        int colsum_sum = 0;
                        int rowsum_conjugated_sum = 0;
                        int k2 = 0;
                        for (; k2 < ncol; k2++) {
                            colsum_sum += seq.colsum[k2];
                            rowsum_conjugated_sum += rowsum_conjugate[k2];
                            if (colsum_sum == rowsum_conjugated_sum) {
                                break;
                            }
                        }

                        // if critical point k2 does exist
                        if (k2 < ncol - 1) {

                            // sequence can be decomposed by splitting column indices at position k2
                            Instance sub1(nrow, k2 + 1);
                            Instance sub2(nrow, ncol - k2 - 1);
                            memcpy(&sub1.colsum[0], &seq.colsum[0], (k2 + 1) * sizeof(int));
                            memcpy(&sub2.colsum[0], &seq.colsum[k2 + 1], (ncol - k2 - 1) * sizeof(int));
                            memcpy(&sub1.colindex[0], &seq.colindex[0], (k2 + 1) * sizeof(int));
                            memcpy(&sub2.colindex[0], &seq.colindex[k2 + 1], (ncol - k2 - 1) * sizeof(int));
                            conjugate(&sub1.rowsum[0], &rowsum_conjugate[0], nrow, k2 + 1);
                            conjugate(&sub2.rowsum[0], &rowsum_conjugate[k2 + 1], nrow, ncol - k2 - 1);
                            memcpy(&sub1.rowindex[0], &seq.rowindex[0], nrow * sizeof(int));
                            memcpy(&sub2.rowindex[0], &seq.rowindex[0], nrow * sizeof(int));

                            decompose_recursive(sub1, f, true, false);
                            decompose_recursive(sub2, f, true, true);
                        }
                        else { // critical point k2 does not exist

                            decompose_recursive(seq, f, false, false);
                        }

                        delete[] rowsum_conjugate;

                    }
                    else { // seq is primitive
                        f(seq); // process sequence
                    }
                }


                /**
                 * (Helper Function) Recursively decompose sequence seq.
                 * For each primitive sequence generated during the process,
                 * the function f is called.
                 * @param seq Realizable sequence pair.
                 * @param f Function called for each primitive sequence.
                 */
                void decompose_recursive(
                        const Instance &seq,
                        const std::function<void( const Instance &

                )> &f) {

                    const int nrow = (int) seq.getNumRows();
                    const int ncol = (int) seq.getNumCols();

                    // std::cout << seq.toString() << std::endl;

                    // build conjugate sequence
                    int *rowsum_conjugate = new int[ncol];
                    int *colsum_conjugate = new int[nrow];
                    conjugate(rowsum_conjugate, &seq.rowsum[0], ncol, nrow);
                    conjugate(colsum_conjugate, &seq.colsum[0], nrow, ncol);

                    // determine critical point k1
                    int rowsum_sum = 0;
                    int colsum_conjugated_sum = 0;
                    int k1 = 0;
                    for (; k1 < nrow; k1++) {
                        rowsum_sum += seq.rowsum[k1];
                        colsum_conjugated_sum += colsum_conjugate[k1];
                        if (rowsum_sum == colsum_conjugated_sum) {
                            break;
                        }
                    }

                    // determine critical point k2
                    int colsum_sum = 0;
                    int rowsum_conjugated_sum = 0;
                    int k2 = 0;
                    for (; k2 < ncol; k2++) {
                        colsum_sum += seq.colsum[k2];
                        rowsum_conjugated_sum += rowsum_conjugate[k2];
                        if (colsum_sum == rowsum_conjugated_sum) {
                            break;
                        }
                    }

                    // if seq is primitive
                    if (k1 == nrow - 1 && k2 == ncol - 1) {
                        f(seq);
                    } else { // seq is composite

                        if (k1 < nrow - 1) {

                            // sequence can be decomposed by splitting row indices at position k1
                            Instance sub1(k1 + 1, ncol);
                            Instance sub2(nrow - k1 - 1, ncol);
                            memcpy(&sub1.rowsum[0], &seq.rowsum[0], (k1 + 1) * sizeof(int));
                            memcpy(&sub2.rowsum[0], &seq.rowsum[k1 + 1], (nrow - k1 - 1) * sizeof(int));
                            memcpy(&sub1.rowindex[0], &seq.rowindex[0], (k1 + 1) * sizeof(int));
                            memcpy(&sub2.rowindex[0], &seq.rowindex[k1 + 1], (nrow - k1 - 1) * sizeof(int));
                            conjugate(&sub1.colsum[0], &colsum_conjugate[0], ncol, k1 + 1);
                            conjugate(&sub2.colsum[0], &colsum_conjugate[k1 + 1], ncol, nrow - k1 - 1);
                            memcpy(&sub1.colindex[0], &seq.colindex[0], ncol * sizeof(int));
                            memcpy(&sub2.colindex[0], &seq.colindex[0], ncol * sizeof(int));

                            decompose_recursive(sub1, f);
                            decompose_recursive(sub2, f);
                        } else { // k2 < ncol - 1

                            // sequence can be decomposed by splitting column indices at position k2
                            Instance sub1(nrow, k2 + 1);
                            Instance sub2(nrow, ncol - k2 - 1);
                            memcpy(&sub1.colsum[0], &seq.colsum[0], (k2 + 1) * sizeof(int));
                            memcpy(&sub2.colsum[0], &seq.colsum[k2 + 1], (ncol - k2 - 1) * sizeof(int));
                            memcpy(&sub1.colindex[0], &seq.colindex[0], (k2 + 1) * sizeof(int));
                            memcpy(&sub2.colindex[0], &seq.colindex[k2 + 1], (ncol - k2 - 1) * sizeof(int));
                            conjugate(&sub1.rowsum[0], &rowsum_conjugate[0], nrow, k2 + 1);
                            conjugate(&sub2.rowsum[0], &rowsum_conjugate[k2 + 1], nrow, ncol - k2 - 1);
                            memcpy(&sub1.rowindex[0], &seq.rowindex[0], nrow * sizeof(int));
                            memcpy(&sub2.rowindex[0], &seq.rowindex[0], nrow * sizeof(int));

                            decompose_recursive(sub1, f);
                            decompose_recursive(sub2, f);

                        }
                    }

                    delete[] rowsum_conjugate;
                    delete[] colsum_conjugate;
                }

            public:


                /**
                 * Create a Decompositor object for a certain instance.
                 * @param inst Row and column sums.
                 */
                Decompositor(const Instance &inst) {

                    _realizable = isRealizable(inst);

                    if (!_realizable)
                        return;

                    const int nrow = (int) inst.getNumRows();
                    const int ncol = (int) inst.getNumCols();

                    // create initial sequence
                    _inst = Instance(nrow, ncol);

                    // sort row and column sums descendingly
                    struct A {
                        int index;
                        int value;
                    };

                    A *rows = new A[nrow];
                    A *cols = new A[ncol];
                    for (int i = 0; i < nrow; i++)
                        rows[i] = {inst.rowindex[i], inst.rowsum[i]};
                    for (int j = 0; j < ncol; j++)
                        cols[j] = {inst.colindex[j], inst.colsum[j]};

                    std::sort(rows, rows + nrow, [](const A &a, const A &b) {
                        return a.value > b.value;
                    });
                    std::sort(cols, cols + ncol, [](const A &a, const A &b) {
                        return a.value > b.value;
                    });

                    for (int i = 0; i < nrow; i++) {
                        _inst.rowindex[i] = rows[i].index;
                        _inst.rowsum[i] = rows[i].value;
                    }

                    for (int j = 0; j < ncol; j++) {
                        _inst.colindex[j] = cols[j].index;
                        _inst.colsum[j] = cols[j].value;
                    }

                    delete[] rows;
                    delete[] cols;
                }


                /**
                  * Decompose the sequence pair into primitives.
                  * The function f is called for each primitive.
                  * @param m Sequence of row and column sums.
                  * @param f Function to be called for each primtive.
                  */
                inline
                void decompose(const std::function<void( const Instance &)> &f) {

                    if (!_realizable)
                        return;

                    // start decomposing the sequence
                    //decompose_recursive(_inst, f);
                    decompose_recursive(_inst, f, true, true);
                }

            };
        }
    }
}

#endif //MARATHON_FIXED_MARGIN_DECOMPOSITION_H