/*
 * quasi_curveball6.h
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

#ifndef QUASI_CURVEBALL6_H_
#define QUASI_CURVEBALL6_H_

#include "marathon/binary_matrix/fixed_margin/markov_chain.h"

namespace marathon {
    namespace binary_matrix {
        namespace fixed_margin {

            /**
             * Markov chain defined in
             *
             * TODO
             * Does the first then the second trading parts
             * Uses THE SAME random bits for the second part
             * Does a second round of trades with a whatever new mask
             */
            class QuasiCurveball6 : public MarkovChain {

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
                    int nTradeBits = nTradableCols / 2;
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

                    return static_cast<int>(std::pow(2, nTradeBits));
                }

                /**
                 * Apply first part of Paired Trades to adjacency matrix M at rows r1 and r2.
                 * @param M
                 * @param r1
                 * @param r2
                 * @param tradeBits
                 */
                int applyFirstTrade(BinaryMatrix &M, const int r1, const int r2, std::bitset<32> tradeBits) const {
                    const int ncol = (int) _inst.getNumCols();

                    const int NOT_FOUND = ncol;

                    int tradeCount = 0;
//                    std::cout << "Trade bits part one: ";

                    int lastTradableCol = NOT_FOUND;
                    for (int j = 0; j < ncol; j++) {

                        int r1_j = static_cast<bool>(M.get(r1, j));
                        int r2_j = static_cast<bool>(M.get(r2, j));

                        if (r1_j != r2_j) {

                            if (lastTradableCol == NOT_FOUND) { lastTradableCol = j; }
                            else {
                                bool doTrade = tradeBits[tradeCount];
//                                std::cout << tradeBits[tradeCount];
                                tradeCount++;
                                if (doTrade) {
                                    bool r1_last = static_cast<bool>(M.get(r1, lastTradableCol));
                                    bool r2_last = static_cast<bool>(M.get(r2, lastTradableCol));

                                    M.set(r1, lastTradableCol, r1_j);
                                    M.set(r2, lastTradableCol, r2_j);

                                    M.set(r1, j, r1_last);
                                    M.set(r2, j, r2_last);
                                }
                                lastTradableCol = NOT_FOUND;
                            }
                        }
                    }
//                    std::cout << "\n";

                    return tradeCount;

                }

                /**
                  * Apply second part of Paired Trades to adjacency matrix M at rows r1 and r2.
                  * @param M
                  * @param r1
                  * @param r2
                  * @param tradeBits
                  */
                int applySecondTrade(BinaryMatrix &M, const int r1, const int r2, std::bitset<32> tradeBits) const {
                    const int ncol = (int) _inst.getNumCols();

                    const int NOT_FOUND = ncol;

                    int tradeCount = 0;
                    int offset = 0; // using the same random bits //getNumTradableCols(M, r1, r2) / 2;

//                    std::cout << "Trade bits part two: ";

                    int lastTradableCol = NOT_FOUND;
                    bool isFirstTradable = true;
                    for (int j = 0; j < ncol; j++) {

                        int r1_j = static_cast<bool>(M.get(r1, j));
                        int r2_j = static_cast<bool>(M.get(r2, j));

                        if (r1_j != r2_j) {
                            if (isFirstTradable) { isFirstTradable = false; }
                            else if (lastTradableCol == NOT_FOUND) { lastTradableCol = j; }
                            else {
                                bool doTrade = tradeBits[offset + tradeCount];
//                                std::cout << tradeBits[offset + tradeCount];
                                tradeCount++;
                                if (doTrade) {
                                    bool r1_last = static_cast<bool>(M.get(r1, lastTradableCol));
                                    bool r2_last = static_cast<bool>(M.get(r2, lastTradableCol));

                                    M.set(r1, lastTradableCol, r1_j);
                                    M.set(r2, lastTradableCol, r2_j);

                                    M.set(r1, j, r1_last);
                                    M.set(r2, j, r2_last);
                                }
                                lastTradableCol = NOT_FOUND;
                            }
                        }
                    }
//                    std::cout << "\n";

                    return tradeCount;

                }

                /**
                  * Apply second part of Paired Trades to adjacency matrix M at rows r1 and r2.
                  * @param M
                  * @param r1
                  * @param r2
                  * @param tradeBits
                  */
                void applyTrade(BinaryMatrix &M, const int r1, const int r2, std::bitset<32> tradeBits) const {
                    const int ncol = (int) _inst.getNumCols();

                    const int NOT_FOUND = ncol;

//                    std::cout << firstOrSecond << "\n";
//                    std::cout << tradeBits << "\n";
//                    std::cout << M.fancyString() << "\n";

                    int tradeCount = 0;
                    tradeCount += applyFirstTrade(M, r1, r2, tradeBits);
                    tradeCount += applySecondTrade(M, r1, r2, tradeBits);

//                    std::cout << M.fancyString() << "\n";

                    // Trades cannot change the number of tradable columns
//                    const int nTradesExpected = std::max(0, getNumTradableCols(M, r1, r2) - 1);
//                    std::cout << "tradeCount " << tradeCount << " nTradesExpected " << nTradesExpected << "\n";
//                    assert(tradeCount == nTradesExpected);
                }

                /**
                  * Returns the single row notation.
                  * @param M
                  * @param r1
                  * @param r2
                 */
                std::bitset<32> getTradableBitsFromA(BinaryMatrix &M, const int r1, const int r2) const {
                    const int ncol = (int) _inst.getNumCols();

                    std::bitset<32> tradableBitsFromA(0);
                    int nTradableCols = 0;
                    for (int j = 0; j < ncol; j++) {

                        int r1_j = static_cast<bool>(M.get(r1, j));
                        int r2_j = static_cast<bool>(M.get(r2, j));

                        if (r1_j != r2_j) {
                            tradableBitsFromA[nTradableCols] = r1_j;
                            nTradableCols++;
                        }
                    }

                    return tradableBitsFromA;
                }



            protected:

                // auxiliary array
                std::vector<int> tmp1;

            public:

                /**
                 * Create a Markov chain.
                 * @param inst Row and column sums.
                 */
                explicit QuasiCurveball6(Instance inst) : MarkovChain(std::move(inst)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                * Create a Markov chain.
                * @param m Binary matrix used as initial state
                */
                explicit QuasiCurveball6(BinaryMatrix m) : MarkovChain(std::move(m)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                 * Create a Markov chain for the given instance.
                 * Use a specified state as initial state.
                 * @param inst Row and Column sums.
                 * @param bin BinaryMatrix used as initial state.
                 */
                QuasiCurveball6(Instance inst, BinaryMatrix bin)
                        : MarkovChain(std::move(inst), std::move(bin)) {
                    tmp1.resize(_currentState.getNumCols());
                }

                /**
                 * Create a Markov chain.
                 * @param inst String-encoded instance.
                 * Instances have the form "2,2,2;1,2,1,2".
                 * The semicolon separates the row sums from the column sums.
                 */
                explicit QuasiCurveball6(const std::string &inst) : QuasiCurveball6(Instance(inst)) {

                }

                /**
                * Create a Markov chain.
                * @param rowsum Sequence of row sums.
                * @param colsum Sequence of column sums.
                */
                QuasiCurveball6(
                        const std::vector<int> &rowsum,
                        const std::vector<int> &colsum
                ) : QuasiCurveball6(Instance(rowsum, colsum)) {

                }

                /**
                 * Create a Markov chain.
                 * @param rowsum Sequence of row sums.
                 * @param colsum Sequence of column sums.
                 * @param nrow Number of rows.
                 * @param ncol Number of columns.
                 */
                QuasiCurveball6(
                        const int *rowsum,
                        const int *colsum,
                        size_t nrow,
                        size_t ncol
                ) : QuasiCurveball6(Instance(rowsum, colsum, nrow, ncol)) {

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
                            const Rational p(1, num_row_sel * nPossibleEndStates);
                            const int nTradeBits = getNumTradeBits(A, i, k);


//                            std::cout << nPossibleEndStates << "\n";
                            for (int trade = 0; trade < nPossibleEndStates; trade++ ) {
                                // simulate 0: U or 1: D trades

                                std::bitset<32> tradeBits(trade);
                                applyTrade(A, i, k, tradeBits);

                                std::cout << A.fancyString() << "\n";
                                auto tradableBitsFromA = getTradableBitsFromA(A, i, k);
                                std::cout << i << " " << k << "\n" << tradableBitsFromA << "\n";

                                tradeBits ^= tradableBitsFromA;
                                applyTrade(A, i, k, ~tradeBits);
//                                tradeBits[nTradeBits] = tradeBits[0];
//                                tradeBits >>= 1;
//                                applyTrade(A, i, k, tradeBits);


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
                    return std::make_unique<QuasiCurveball6>(_inst, _currentState);
                }
            };
        }
    }
}


#endif /* QUASI_CURVEBAll6_H_ */
