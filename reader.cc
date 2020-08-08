#include "reader.hh"

#include <iostream>

#include "utils.hh"

// Suppress optimization to speed up compilation
#if defined(__GNUC__) && __GNUC__ >=4
#pragma GCC optimize "-O"
#endif

// Suppress spurious warning triggered by toml11 for gcc 4.8.5
#if defined(__GNUC__) && __GNUC__ < 5
#pragma GCC ignored "-Wno-unused-parameter"
#endif

#include "toml11/toml.hpp"

// Re-enable warnings
#if defined(__GNUC__) && __GNUC__ < 5
#pragma GCC warning "-Wno-unused-parameter"
#endif

class Node {
public:
	typedef toml::basic_value<toml::discard_comments,
														std::unordered_map, std::vector> value_type;
	Node(const value_type& v) : v(v) {}
	value_type v;
};

using namespace std;

InputNode::InputNode(const Node n, const std::string& label) :
		n_(std::shared_ptr<Node>(new Node{n})), label_(label) {}

size_t InputNode::size() const
{
	return n().v.size();
}
const Node& InputNode::n() const
{
	return *n_;
}
bool InputNode::contains(const std::string& name) const
{
	return n().v.contains(name);
}
// Inherit label from parent -- could potentially add index here
InputNode InputNode::operator[](
	unsigned int index) const
{
	return InputNode(Node(toml::find(n().v,index)),label());
}
// Could potentially include parent information here
InputNode InputNode::operator[](
	const std::string& s) const
{
	return InputNode(Node(toml::find(n().v,s)),s);
}
std::string InputNode::stringfield_unchecked(
	const std::string& name) const
{
	return toml::find<std::string>(n().v,name);
}
std::string InputNode::stringfield(
	const char *name) const
{
	if(!contains(name)) {
		ostringstream oss;
		oss << "A '" << name <<
			"' property must be specified in '" << label() << "'.";
		emsgroot(oss.str().c_str());
	}
	return stringfield_unchecked(name);
}
double InputNode::numberfield_unchecked(
	const std::string& name) const
{
	const auto value = toml::find(n().v,name);
	if(value.is_floating())
		return value.as_floating();
	else
		return value.as_integer();
}
double InputNode::numberfield(
	const char *title,
	const char *name) const
{
	if(!contains(name)) {
		ostringstream oss;
		oss << "A numeric value for '" << name <<
			"' must be specified in '" << title << "'.";
		emsgroot(oss.str().c_str());
	}
	return numberfield_unchecked(name);
}
int InputNode::intfield_unchecked(
	const std::string& name) const
{
	return toml::find<int>(n().v,name);
}
/// Gets a list of all the keys
vector<string> InputNode::get_keys( ) const
{
	vector<string> keys;
	for(const auto& p : n().v.as_table())
	{
		keys.push_back(p.first);
	}
	return keys;
}

/// Factory to generate InputNode from names file
InputNode parsefile(const std::string& inputfilename)
{
	Node n{toml::parse(inputfilename)};
	// Allow using values from another TOML file as a base for this one. TODO:
	// make this into functions so you can do this recursively.
	if (n.v.contains("baseinputfile")) {
		const string basefile = toml::find<std::string>(n.v,"baseinputfile");
		// TODO: make the filename relative to the original TOML file
		decltype(toml::parse(basefile)) basetomlddata = toml::parse(basefile);
		
		for(const auto& p : basetomlddata.as_table())
		{
			if (!n.v.contains(p.first)) {
				n.v[p.first] = p.second;
			}
		}
	}
	return InputNode(n,"");
}

InputData::InputData(const std::string& inputfilename) :
	data(parsefile(inputfilename))
{}

/// Gets a list of all the keys
vector<string> InputData::get_keys( ) const
{
	return data.get_keys();
}
