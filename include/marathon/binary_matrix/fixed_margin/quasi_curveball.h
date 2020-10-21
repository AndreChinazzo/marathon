/*
 * quasi_curveball.h
 *
 * Created on: Oct 13, 2019
 * Author: Andre Lucas Chinazzo <chinazzo@eit.uni-kl.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
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

#ifndef QUASI_CURVEBALL_H_
#define QUASI_CURVEBALL_H_

#include "marathon/binary_matrix/fixed_margin/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Markov chain defined in
             *
             * TODO
             */
            class QuasiCurveball : public MarkovChain {

            private:

                /**
                 * Returns the number of tradable columns between r1 and r2.
                 * @param M
                 * @param r1
                 * @param r2
                 */
                int getNumTradableCols(BinaryMatrix &M, const int r1, const int r2) const {
                    const int ncol = (int) _inst.getNumCols();

                    int nTradableCols = 0;
                    for (int j = 0; j < ncol; j++) {

                        int r1_j = static_cast<bool>(M.get(r1, j));
                        int r2_j = static_cast<bool>(M.get(r2, j));

                        if (r1_j != r2_j) {
                            nTradableCols++;
                        }
                    }

                    return nTradableCols;
                }

                /**
                 * Returns the number of col-pair-wise trades decisions between r1 and r2.
                 * @param M
                 * @param r1
                 * @param r2
                 */
                int getNumTradeBits(BinaryMatrix &M, const int r1, const int r2) const {

                    int nTradableCols = getNumTradableCols(M, r1, r2);

                    int nTradeBits = std::max(0, nTradableCols-1);
                    return nTradeBits;
                }

                /**
                 * Returns the number of outgoing edges if r1 and r2 are chosen.
                 * @param M
                 * @param r1
                 * @param r2
                 */
                int getNumEndStates(BinaryMatrix &M, const int r1, const int r2) const {

                    int nTradeBits = getNumTradeBits(M, r1, r2);

                    return std::max(1, static_cast<int>(std::pow(2, nTradeBits)));
                }

                /**
                 * Apply Paired Trades to adjacency matrix M at rows r1 and r2.
                 * @param M
                 * @param r1
                 * @param r2
                 * @param tradeBits
                 */
                void applyTrade(BinaryMatrix &M, const int r1, const int r2, std::bitset<32> tradeBits) const {
                    const int ncol = (int) _inst.getNumCols();

                    const int NOT_FOUND = ncol;

                    int tradeCount = 0;

//                    std::cout << tradeBits << "\n";
//                    std::cout << M.fancyString() << "\n";

                    // First pass
                    int lastTradableCol = NOT_FOUND;
                    for (int j = 0; j < ncol; j++) {

                        int r1_j = static_cast<bool>(M.get(r1, j));
                        int r2_j = static_cast<bool>(M.get(r2, j));

                        if (r1_j != r2_j) {

                            if ( lastTradableCol == NOT_FOUND ) { lastTradableCol = j; }
                            else if ( tradeBits[tradeCount++] ) {
                                bool r1_last = static_cast<bool>(M.get(r1, lastTradableCol));
                                bool r2_last = static_cast<bool>(M.get(r2, lastTradableCol));

                                M.set(r1, lastTradableCol, r1_j);
                                M.set(r2, lastTradableCol, r2_j);

                                M.set(r1, j, r1_last);
                                M.set(r2, j, r2_last);

                                lastTradableCol = NOT_FOUND;
                            }
                            else { lastTradableCol = NOT_FOUND; }
                        }
                    }

//                    std::cout << M.fancyString() << "\n";

                    // Second pass, notice the j=1
                    lastTradableCol = NOT_FOUND;
                    bool isFirstTradable = true;
                    for (int j = 0; j < ncol; j++) {

                        int r1_j = static_cast<bool>(M.get(r1, j));
                        int r2_j = static_cast<bool>(M.get(r2, j));

                        if (r1_j != r2_j) {
                            if ( isFirstTradable ) { isFirstTradable = false; }
                            else if ( lastTradableCol == NOT_FOUND ) { lastTradableCol = j; }
                            else if ( tradeBits[tradeCount++] ) {
                                bool r1_last = static_cast<bool>(M.get(r1, lastTradableCol));
                                bool r2_last = static_cast<bool>(M.get(r2, lastTradableCol));

                                M.set(r1, lastTradableCol, r1_j);
                                M.set(r2, lastTradableCol, r2_j);

                                M.set(r1, j, r1_last);
                                M.set(r2, j, r2_last);

                                lastTradableCol = NOT_FOUND;
                            }
                            else { lastTradableCol = NOT_FOUND; }
                        }
                    }

//                    std::cout << M.fancyString() << "\n";

                    // Used them all
//                    std::cout << "tradeCount" << tradeCount << "getNumTradeBits" << getNumTradeBits(M, r1, r2) << "\n";
                    assert(tradeCount == getNumTradeBits(M, r1, r2));
                }

            protected:

                // auxiliary array
                std::vector<int> tmp1;

            public:

                /**
                 * Create a Markov chain.
                 * @param inst Row and column sums.
                 */
                explicit QuasiCurveball(Instance inst) : MarkovChain(std::move(inst)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                * Create a Markov chain.
                * @param m Binary matrix used as initial state
                */
                explicit QuasiCurveball(BinaryMatrix m) : MarkovChain(std::move(m)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                 * Create a Markov chain for the given instance.
                 * Use a specified state as initial state.
                 * @param inst Row and Column sums.
                 * @param bin BinaryMatrix used as initial state.
                 */
                QuasiCurveball(Instance inst, BinaryMatrix bin)
                        : MarkovChain(std::move(inst), std::move(bin)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                 * Create a Markov chain.
                 * @param inst String-encoded instance.
                 * Instances have the form "2,2,2;1,2,1,2".
                 * The semicolon separates the row sums from the column sums.
                 */
                explicit QuasiCurveball(const std::string &inst) : QuasiCurveball(Instance(inst)) {

                }

                /**
                * Create a Markov chain.
                * @param rowsum Sequence of row sums.
                * @param colsum Sequence of column sums.
                */
                QuasiCurveball(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : QuasiCurveball(Instance(rowsum, colsum)) {

                }

                /**
                 * Create a Markov chain.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
                 */
                QuasiCurveball(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) : QuasiCurveball(Instance(rowsum, colsum, nrow, ncol)) {

                }

                /**
                 * Generate each adjacent state x to s and the corresponding proposal propability p(s,x).
                 * For each pair (x,p) call the function f.
                 * @param s
                 * @param process
                 */
                virtual void adjacentStates(
                        const State &x,
                        const std::function<void(const State &, const marathon::Rational &)> &f
                ) const override {

                    const int nrow = (int) _inst.getNumRows();
                    const int ncol = (int) _inst.getNumCols();

                    // convert state reference
                    const BinaryMatrix &X = static_cast<const BinaryMatrix &>(x);

                    // Make a copy of the current state
                    BinaryMatrix A(X);

                    /**
                     * TODO: UPDATE
                     * Definition
                     */

                    // randomly select two row indices
                    for (int i = 0; i < nrow; i++) {
                        for (int k = i + 1; k < nrow; k++) {

                            // calculate the probability of this choice
                            const Integer num_row_sel = binom(nrow, 2);
                            const Integer nPossibleEndStates = getNumEndStates(A, i, k);
                            const Rational p(1, num_row_sel * nPossibleEndStates); // 2: either D or U trade

//                            std::cout << nPossibleEndStates << "\n";
                            for (int trade = 0; trade < nPossibleEndStates; trade++ ) {
                                // simulate 0: U or 1: D trades

                                std::bitset<32> tradeBits(trade);
                                applyTrade(A, i, k, tradeBits);
                                assert(_inst.isValid(A));

                                // process adjacent state
                                f(A, p);

                                // restore original state of rows i and k
                                for (int j = 0; j < ncol; j++) {
                                    A.set(i, j, X.get(i, j));
                                    A.set(k, j, X.get(k, j));
                                }
                            }
                        }
                    }
                }

                /**
                 * Randomize the current state of the Markov chain.
                 */
                virtual void step() override {

                    // TODO: UNTESTED
                    const size_t nrow = _inst.getNumRows();
                    const size_t ncol = _inst.getNumCols();

                    if (nrow == 1 || ncol == 1)
                        return;

                    // randomly select two row indices
                    int i = rg.nextInt(nrow);
                    int k = rg.nextInt(nrow);
                    while (i == k)
                        k = rg.nextInt(nrow);

                    int nTradeBits = getNumTradeBits(_currentState, i, k);

                    std::bitset<32> tradeBits(rg.nextInt(nTradeBits));

                    applyTrade(_currentState, i, k, tradeBits);
                    assert(_inst.isValid(_currentState));
                }

                /**
                 * Create a copy of this MarkovChain.
                 * @return
                 */
                virtual std::unique_ptr<marathon::MarkovChain> copy() const override {
                    return std::make_unique<QuasiCurveball>(_inst, _currentState);
                }
            };
        }
    }
}


#endif /* QUASI_CURVEBALL_H_ */
