#pragma once


struct Locus
{
	virtual ui get_start() = 0;
	virtual string to_string_debug() = 0;
	virtual string to_string_summary() = 0;
	virtual bool is_crispr() = 0;
	virtual bool is_domain() = 0;
	virtual bool is_gene() = 0;
};


