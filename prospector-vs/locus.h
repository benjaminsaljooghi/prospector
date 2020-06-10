#pragma once


struct Locus
{
	virtual ui get_start() = 0;
	virtual string to_string_debug() = 0;
	virtual string to_string_summary() = 0;
};