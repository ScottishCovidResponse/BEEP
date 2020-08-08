#ifndef BEEPMBP__READER_HH
#define BEEPMBP__READER_HH

#include <vector>
#include <memory>

class Node;

class InputNode {
	std::shared_ptr<Node> n_;
	std::string label_;
public:
	explicit InputNode(const Node n, const std::string& label);
	size_t size() const;
	const std::string& label() const
		{
			return label_;
		}
	bool contains(const std::string& name) const;
	InputNode operator[](unsigned int index) const;
	InputNode operator[](const std::string& s) const;
	std::string stringfield_unchecked(
		const std::string& name) const;
	std::string stringfield(
		const char *name) const;
	double numberfield_unchecked(
		const std::string& name) const;
	double numberfield(
		const char *title,
		const char *name) const;
	int intfield_unchecked(
		const std::string& name) const;
	std::vector<std::string> keys() const;
private:
	const Node& n() const;
};

/// Factory to generate InputNode from names file
InputNode parsefile(const std::string& inputfilename);

class InputData {
public:
	InputData(const std::string& inputfilename) :
		data(parsefile(inputfilename)) {}
	InputData(const InputData& data) = delete;
	InputData(InputData&& data) = delete;
	InputData& operator=(const InputData& data) = delete;
	InputData& operator=(InputData&& data) = delete;
	~InputData() {}
	bool contains(const std::string& name) const
		{
			return data.contains(name);
		}
	InputNode open(const std::string& name)
		{
			return data[name];
		}
	/// Gets a list of all the keys
	std::vector<std::string> keys() const
		{
			return data.keys();
		}
	// Information from the TOML file	
	InputNode data;
};

#endif
