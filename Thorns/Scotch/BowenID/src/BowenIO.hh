#ifndef BOWEN_IO_HEADER
#define BOWEN_IO_HEADER
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

/*
class OutputTable
{
	private:
		int width;
		std::ostream& stream_obj;

	public:
		OutputTable(std::ostream& obj, int w): width(w), stream_obj(obj) {}

		template<typename T>
			OutputTable& operator<<(const T& output)
			{
				stream_obj << std::setw(width) << output;

				return *this;
			}

		OutputTable& operator<<(std::ostream& (*func)(std::ostream&))
		{
			func(stream_obj);
			return *this;
		}
};
*/

class OutputTable
{
	private:
		int width;
		bool testsuite;
		std::ofstream stream_obj;

	public:
		OutputTable(std::stringstream &fname, int w, bool tsuite);

		~OutputTable();

		template<typename T>
			OutputTable& operator<<(const T& output);

		OutputTable& operator<<(std::ostream& (*func)(std::ostream&));
};

OutputTable::OutputTable(std::stringstream &fname, int w, bool tsuite)
{
	stream_obj.open(fname.str().c_str());
	width = w;
	testsuite = tsuite;
}

OutputTable::~OutputTable()
{
	stream_obj.close();
}


template<typename T>
OutputTable::OutputTable& OutputTable::operator<<(const T& output)
{
	if(testsuite)
		stream_obj << std::setw(width) << output;
	return *this;	
}

OutputTable::OutputTable& OutputTable::operator<<(std::ostream& (*func)(std::ostream&))
{
	if(testsuite)
		func(stream_obj);
	return *this;
}
#endif

