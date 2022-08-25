#pragma once

#ifndef PROSPECTOR_CAS_HELPERS_H
#define PROSPECTOR_CAS_HELPERS_H

#include "cas.h"

static unordered_set<string> alternate_starts_pos { "GTG", "TTG" };
static unordered_set<string> alternate_starts_neg { "CAC", "CAA" };

static unordered_set<string> start_codons_pos{ "ATG" };
static unordered_set<string> start_codons_neg{ "CAT" };
static unordered_set<string> stop_codons_pos{ "TAA", "TAG", "TGA" };
static unordered_set<string> stop_codons_neg{ "TTA", "CTA", "TCA" };
static ull steps = 200;

ull generic_expand(const string& genome, unordered_set<string>& sinks, ull begin, ull increment, ull increments);

void fragment_expansion_pos(Fragment* fragment, const string& genome);

void fragment_expansion_neg(Fragment* fragment, const string& genome);

vector<vector<ull>> cluster_index(const vector<ull>& indices);

bool good_clusters(const vector<vector<ull>>& clusters);

ull demarc_start_clusters(const vector<vector<ull>>& clusters);

ull demarc_final_clusters(const vector<vector<ull>>& clusters);

void compute_demarc(Fragment* frag);

void compute_details(Fragment* fragment, const string& genome);

bool fragment_equivalent(const Fragment* a, const Fragment* b);

bool fragment_contains(const Fragment* a, const Fragment* b);

::size_t allocate_profiles_on_gpu(
        vector<CasProfile *> &profiles,
        ui *&h_profile_kmers, ui *&d_profile_kmers,
        ui *&h_profile_coords, ui *&d_profile_coords);

::size_t allocate_translations_on_gpu(
        vector<Translation *> &translations,
        kmer *&h_translation_kmers, kmer *&d_translation_kmers,
        ui *&h_translation_coords, ui *&d_translation_coords);

#endif //PROSPECTOR_CAS_HELPERS_H
