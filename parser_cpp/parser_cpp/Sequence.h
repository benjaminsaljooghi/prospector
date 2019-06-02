#pragma once

#include <vector>
#include <string>
#include "consts.h"

class Sequence
{
    std::string seq;
    int start_pos;

    public:

    Sequence(std::string seq, int start);

    CUDA_CALLABLE_MEMBER int length();

    CUDA_CALLABLE_MEMBER int start();

    CUDA_CALLABLE_MEMBER int end();

    CUDA_CALLABLE_MEMBER Sequence subseq(int start, int length);

    CUDA_CALLABLE_MEMBER std::string sequence();

    CUDA_CALLABLE_MEMBER char operator[](int i);

    CUDA_CALLABLE_MEMBER bool is_dyad();

    CUDA_CALLABLE_MEMBER std::vector<Sequence> dyads(int k);

    CUDA_CALLABLE_MEMBER std::vector<Sequence> dyads(int k_start, int k_end);

    CUDA_CALLABLE_MEMBER bool operator<(const Sequence rhs) const;
};

