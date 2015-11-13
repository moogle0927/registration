#ifndef __INFO__
#define __INFO__

#include "common.h"
#include "nariinfocontroller.h"
#include "narifile.h"
#include "naricommon.h"
#include <string>

struct text_info
{
	
	std::string dirLabel;         
	std::string set;          
	std::string pathId;          
	std::string case_list; 
	std::string dirOrg ;
	std::string dirDisp;
	std::string dirO ;
	std::string dirBase ;
	nari::case_t caseRef;
	int tmp ;
	int rangex ;
	int rangey ;
	int rangez ;
	int d ;

					
	inline void input(const std::string &path) // テキストから入力情報を書き込み
	{
		nari::infocontroller info;
		info.load(path);

		dirLabel = nari::file::add_delim(info.get_as_str("in"));
		set = info.get_as_str("set");
		pathId = nari::file::add_delim(info.get_as_str("case"));
		pathId += set;
		dirOrg = nari::file::add_delim(info.get_as_str("out"))+info.get_as_str("01-01");
		dirDisp = nari::file::add_delim(info.get_as_str("out"))+info.get_as_str("07-01");
		dirO = nari::file::add_delim(info.get_as_str("forLearn"))+info.get_as_str("labelStd");
		dirBase = dirDisp;
		caseRef = (info.get_as_str("ref_file"),info.get_as_str("ref_name"));
		
		tmp = info.get_as_int("tmp_size");
		rangex = info.get_as_int("range_x");
		rangey = info.get_as_int("range_y");
		rangez = info.get_as_int("range_z");
		d = info.get_as_int("search_interval");
		
		

		info.output();
	}
};


#endif
