#include "reader.hh"

#include <iostream>

#include "utils.hh"

using namespace std;

size_t InputNode::size() const
{
	return n.size();
}
bool InputNode::contains(const std::string& name) const
{
	return n.contains(name);
}
// Inherit label from parent -- could potentially add index here
InputNode InputNode::operator[](
	unsigned int index) const
{
	return InputNode(toml::find(n,index),label());
}
// Could potentially include parent information here
InputNode InputNode::operator[](
	const std::string& s) const
{
	return InputNode(toml::find(n,s),s);
}
std::string InputNode::stringfield_unchecked(
	const std::string& name) const
{
	return toml::find<std::string>(n,name);
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
	const auto value = toml::find(n,name);
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
	return toml::find<int>(n,name);
}

InputData::InputData(const std::string& inputfilename) :
	data(toml::parse(inputfilename),"")
{
	// Allow using values from another TOML file as a base for this one. TODO:
	// make this into functions so you can do this recursively.
	if (contains("baseinputfile")) {
		const string basefile = data.stringfield_unchecked("baseinputfile");
		// TODO: make the filename relative to the original TOML file
		decltype(toml::parse(basefile)) basetomlddata = toml::parse(basefile);
		
		for(const auto& p : basetomlddata.as_table())
		{
			if (!contains(p.first)) {
				data.n[p.first] = p.second;
			}
		}
	}
}
/// Gets a list of all the keys
vector<string> InputData::get_keys( ) const
{
	vector<string> keys;
	for(const auto& p : data.n.as_table())
	{
		keys.push_back(p.first);
	}
	return keys;
}
