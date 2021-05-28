// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_ATTRIBUTE_H
#define HEPMC3_ATTRIBUTE_H
/**
 *  @file Attribute.h
 *  @brief Definition of \b class Attribute, \b class IntAttribute and \b class StringAttribute
 *
 *  @class HepMC3::Attribute
 *  @brief Base class for all attributes
 *
 *  Contains virtual functions to_string and from_string that
 *  each attribute must implement, as well as init function that
 *  attributes should overload to initialize parsed attribute
 *
 *  @ingroup attributes
 *
 */
#include <cstdio> // sprintf
#include <string>
#include <limits>
#include <sstream>
#include <iomanip>
#include <map>

#include "HepMC3/GenParticle_fwd.h"
#include "HepMC3/GenVertex_fwd.h"

/** Deprecated */
using std::string;

namespace HepMC3 {

/** @brief Forward declaration of GenEvent. */
class GenEvent;

/** @brief Forward declaration of GenRunInfo. */
class GenRunInfo;

/** @brief Forward declaration of GenParticle. */
// class GenParticle;
class Attribute {
//
// Constructors
//
public:
    /** @brief Default constructor */
    //Note: m_event should be set to nullptr in case event is deleted!
    Attribute():m_is_parsed(true) {    m_event=nullptr; }

    /** @brief Virtual destructor */
    virtual ~Attribute() {}

protected:
    /** @brief Protected constructor that allows to set string
     *
     *  Used when parsing attributes from file. An StringAttribute class
     *  object is made, which uses this constructor to signify that
     *  it just holds string without parsing it.
     *
     *  @note There should be no need for user class to ever use this constructor
     */
    //Note: m_event should be set to nullptr n case event is deleted!
    explicit Attribute(const std::string &st):m_is_parsed(false),m_string(st) { m_event=nullptr; }

    /** @brief GenEvent is a friend */
    friend class GenEvent;

//
// Virtual Functions
//
public:
    /** @brief Fill class content from string.
     */
    virtual bool from_string(const std::string & att) = 0;

    /** @brief Optionally initialize the attribute after from_string.
     */
    virtual bool init() {
        return true;
    }

    /** @brief Optionally initialize the attribute after from_string
     *
     * Is passed a reference to the GenRunInfo object to which the
     * Attribute belongs.
     */
    virtual bool init(const GenRunInfo & ) {
        return true;
    }

    /** @brief Fill string from class content */
    virtual bool to_string(std::string &att) const = 0;

//
// Accessors
//
public:
    /** @brief Check if this attribute is parsed */
    bool is_parsed() const { return m_is_parsed; }

    /** @brief Get unparsed string */
    const std::string& unparsed_string() const { return m_string; }

    /** return the GenEvent to which this Attribute belongs, if at all. */
    const GenEvent * event() const {
        return m_event;
    }

    /** return the GenParticle to which this Attribute belongs, if at all. */
    GenParticlePtr particle() {
        return m_particle;
    }

    /** return the GenParticle to which this Attribute belongs, if at all. */
    ConstGenParticlePtr particle() const {
        return std::const_pointer_cast<GenParticle>(m_particle);
    }

    /** return the GenVertex to which this Attribute belongs, if at all. */
    GenVertexPtr vertex() {
        return m_vertex;
    }

    /** return the GenVertex to which this Attribute belongs, if at all. */
    ConstGenVertexPtr vertex() const {
        return std::const_pointer_cast<GenVertex>(m_vertex);
    }

protected:
    /** @brief Set is_parsed flag */
    void set_is_parsed(bool flag) { m_is_parsed = flag; }

    /** @brief Set unparsed string */
    void set_unparsed_string(const std::string &st) { m_string = st; }

//
// Fields
//
private:
    bool   m_is_parsed;             //!< Is this attribute parsed?
    std::string m_string;                //!< Raw (unparsed) string
    const GenEvent * m_event;       //!< Possibility to be aware of the
    //!  controlling GenEvent object.
    GenParticlePtr m_particle; //!< Particle to which assigned.
    GenVertexPtr m_vertex;      //!< Vertex to which assigned.
};

/**
 *  @class HepMC3::IntAttribute
 *  @brief Attribute that holds an Integer implemented as an int
 *
 *  @ingroup attributes
 */
class IntAttribute : public Attribute {
public:

    /** @brief Default constructor */
    IntAttribute():Attribute(),m_val(0) {}

    /** @brief Constructor initializing attribute value */
    IntAttribute(int val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        m_val = atoi( att.c_str() );
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att = std::to_string(m_val);
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    int value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const int& i) {
        m_val = i;
    }

private:
    int m_val; ///< Attribute value
};

/**
 *  @class HepMC3::LongAttribute
 *  @brief Attribute that holds an Integer implemented as an int
 *
 *  @ingroup attributes
 */
class LongAttribute : public Attribute {
public:

    /** @brief Default constructor */
    LongAttribute(): Attribute(), m_val(0) {}

    /** @brief Constructor initializing attribute value */
    LongAttribute(long val): Attribute(), m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        m_val = atol( att.c_str() );
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att = std::to_string(m_val);
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    long value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const long& l) {
        m_val = l;
    }

private:

    long m_val; ///< Attribute value

};

/**
 *  @class HepMC3::DoubleAttribute
 *  @brief Attribute that holds a real number as a double.
 *
 *  @ingroup attributes
 */
class DoubleAttribute : public Attribute {
public:

    /** @brief Default constructor */
    DoubleAttribute(): Attribute(), m_val(0.0) {}

    /** @brief Constructor initializing attribute value */
    DoubleAttribute(double val): Attribute(), m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        m_val = atof( att.c_str() );
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        std::ostringstream oss;
        oss << std::setprecision(std::numeric_limits<double>::digits10)
            << m_val;
        att = oss.str();
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    double value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const double& d) {
        m_val = d;
    }

private:

    double m_val; ///< Attribute value
};

/**
 *  @class HepMC3::FloatAttribute
 *  @brief Attribute that holds a real number as a float.
 *
 *  @ingroup attributes
 */
class FloatAttribute : public Attribute {
public:

    /** @brief Default constructor */
    FloatAttribute(): Attribute(), m_val(0.0) {}

    /** @brief Constructor initializing attribute value */
    FloatAttribute(float val): Attribute(), m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        m_val = float(atof( att.c_str() ));
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        std::ostringstream oss;
        oss << std::setprecision(std::numeric_limits<float>::digits10)
            << m_val;
        att = oss.str();
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    float value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const float& f) {
        m_val = f;
    }

private:

    float m_val; ///< Attribute value
};

/**
 *  @class HepMC3::StringAttribute
 *  @brief Attribute that holds a string
 *
 *  Default attribute constructed when reading input files.
 *  It can be then parsed by other attributes or left as a string.
 *
 *  @ingroup attributes
 *
 */
class StringAttribute : public Attribute {
public:

    /** @brief Default constructor - empty string */
    StringAttribute():Attribute() {}

    /** @brief String-based constructor
     *
     *  The Attribute constructor used here marks that this is an unparsed
     *  string that can be (but does not have to be) parsed
     *
     */
    StringAttribute(const std::string &st):Attribute(st) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        set_unparsed_string(att);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att = unparsed_string();
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::string value() const {
        return unparsed_string();
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::string& s) {
        set_unparsed_string(s);
    }

};

/**
 *  @class HepMC3::CharAttribute
 *  @brief Attribute that holds an Chareger implemented as an int
 *
 *  @ingroup attributes
 */
class CharAttribute : public Attribute {
public:

    /** @brief Default constructor */
    CharAttribute():Attribute(),m_val(0) {}

    /** @brief Constructor initializing attribute value */
    CharAttribute(char val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        if (att.size())
        {
            m_val = att.at(0);
            return true;
        }
        return false;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const override {
        att = std::to_string(m_val);
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    char value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const char& i) {
        m_val = i;
    }

private:
    char m_val; ///< Attribute value
};

/**
 *  @class HepMC3::LongLongAttribute
 *  @brief Attribute that holds an Integer implemented as an int
 *
 *  @ingroup attributes
 */
class LongLongAttribute : public Attribute {
public:

    /** @brief Default constructor */
    LongLongAttribute(): Attribute(), m_val(0) {}

    /** @brief Constructor initializing attribute value */
    LongLongAttribute(long long val): Attribute(), m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        m_val = atoll( att.c_str() );
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att = std::to_string(m_val);
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    long long value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const long long& l) {
        m_val = l;
    }

private:

    long  long m_val; ///< Attribute value

};

/**
 *  @class HepMC3::LongDoubleAttribute
 *  @brief Attribute that holds a real number as a double.
 *
 *  @ingroup attributes
 */
class LongDoubleAttribute : public Attribute {
public:

    /** @brief Default constructor */
    LongDoubleAttribute(): Attribute(), m_val(0.0) {}

    /** @brief Constructor initializing attribute value */
    LongDoubleAttribute(long double val): Attribute(), m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        m_val = strtold( att.c_str(),NULL);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        std::ostringstream oss;
        oss << std::setprecision(std::numeric_limits<long double>::digits10)
            << m_val;
        att = oss.str();
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    long double value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const long double& d) {
        m_val = d;
    }

private:

    long double m_val; ///< Attribute value
};



/**
 *  @class HepMC3::UIntAttribute
 *  @brief Attribute that holds an unsigned int
 *
 *  @ingroup attributes
 */
class UIntAttribute : public Attribute {
public:

    /** @brief Default constructor */
    UIntAttribute():Attribute(),m_val(0) {}

    /** @brief Constructor initializing attribute value */
    UIntAttribute(unsigned int val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        m_val = strtoul(att.c_str(), NULL, 0);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att = std::to_string(m_val);
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    unsigned int value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const unsigned int& i) {
        m_val = i;
    }

private:
    unsigned int m_val; ///< Attribute value
};



/**
 *  @class HepMC3::ULongAttribute
 *  @brief Attribute that holds an unsigned long
 *
 *  @ingroup attributes
 */
class ULongAttribute : public Attribute {
public:

    /** @brief Default constructor */
    ULongAttribute():Attribute(),m_val(0) {}

    /** @brief Constructor initializing attribute value */
    ULongAttribute(unsigned long val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        m_val = strtoul(att.c_str(), NULL, 0);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att = std::to_string(m_val);
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    unsigned long value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const unsigned long& i) {
        m_val = i;
    }

private:
    unsigned long m_val; ///< Attribute value
};


/**
 *  @class HepMC3::ULongLongAttribute
 *  @brief Attribute that holds an unsigned long long
 *
 *  @ingroup attributes
 */
class ULongLongAttribute : public Attribute {
public:

    /** @brief Default constructor */
    ULongLongAttribute():Attribute(),m_val(0) {}

    /** @brief Constructor initializing attribute value */
    ULongLongAttribute(unsigned long long val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        m_val = strtoull(att.c_str(), NULL, 0);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att = std::to_string(m_val);
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    unsigned long long value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const unsigned long long& i) {
        m_val = i;
    }

private:
    unsigned long long m_val; ///< Attribute value
};
/**
 *  @class HepMC3::BoolAttribute
 *  @brief Attribute that holds an Booleger implemented as an int
 *
 *  @ingroup attributes
 */
class BoolAttribute : public Attribute {
public:

    /** @brief Default constructor */
    BoolAttribute():Attribute(),m_val(false) {}

    /** @brief Constructor initializing attribute value */
    BoolAttribute(bool val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override{
        if (att.size()!=1) return false;
        if(att==std::string("1")) {m_val = true;  return true;}
        if(att==std::string("0")) {m_val = false; return true;}
        return false;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const override{
        att = std::to_string(m_val);
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    bool value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const bool& i) {
        m_val = i;
    }

private:
    bool m_val; ///< Attribute value
};

/**
 *  @class HepMC3::VectorCharAttribute
 *  @brief Attribute that holds a vector of charegers of type  char
 *
 *  @ingroup attributes
 */
class VectorCharAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorCharAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorCharAttribute(std::vector<char> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        char  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<char> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<char>& i) {
        m_val = i;
    }

private:
    std::vector<char> m_val; ///< Attribute value
};

/**
 *  @class HepMC3::VectorFloatAttribute
 *  @brief Attribute that holds a vector of floategers of type  float
 *
 *  @ingroup attributes
 */
class VectorFloatAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorFloatAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorFloatAttribute(std::vector<float> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        float  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<float> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<float>& i) {
        m_val = i;
    }

private:
    std::vector<float> m_val; ///< Attribute value
};


/**
 *  @class HepMC3::VectorLongDoubleAttribute
 *  @brief Attribute that holds a vector of long doubleegers of type  long double
 *
 *  @ingroup attributes
 */
class VectorLongDoubleAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorLongDoubleAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorLongDoubleAttribute(std::vector<long double> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        long double  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<long double> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<long double>& i) {
        m_val = i;
    }

private:
    std::vector<long double> m_val; ///< Attribute value
};



/**
 *  @class HepMC3::VectorLongLongAttribute
 *  @brief Attribute that holds a vector of long longegers of type  long long
 *
 *  @ingroup attributes
 */
class VectorLongLongAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorLongLongAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorLongLongAttribute(std::vector<long long> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        long long  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<long long> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<long long>& i) {
        m_val = i;
    }

private:
    std::vector<long long> m_val; ///< Attribute value
};

/**
 *  @class HepMC3::VectorUIntAttribute
 *  @brief Attribute that holds a vector of unsigned integers of type  unsigned int
 *
 *  @ingroup attributes
 */
class VectorUIntAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorUIntAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorUIntAttribute(std::vector<unsigned int> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        unsigned int  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<unsigned int> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<unsigned int>& i) {
        m_val = i;
    }

private:
    std::vector<unsigned int> m_val; ///< Attribute value
};

/**
 *  @class HepMC3::VectorULongAttribute
 *  @brief Attribute that holds a vector of unsigned longegers of type  unsigned long
 *
 *  @ingroup attributes
 */
class VectorULongAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorULongAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorULongAttribute(std::vector<unsigned long> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        unsigned long  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<unsigned long> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<unsigned long>& i) {
        m_val = i;
    }

private:
    std::vector<unsigned long> m_val; ///< Attribute value
};


/**
 *  @class HepMC3::VectorULongLongAttribute
 *  @brief Attribute that holds a vector of unsigned long longegers of type  unsigned long long
 *
 *  @ingroup attributes
 */
class VectorULongLongAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorULongLongAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorULongLongAttribute(std::vector<unsigned long long> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        unsigned long long  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<unsigned long long> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<unsigned long long>& i) {
        m_val = i;
    }

private:
    std::vector<unsigned long long> m_val; ///< Attribute value
};

/**
 *  @class HepMC3::VectorIntAttribute
 *  @brief Attribute that holds a vector of integers of type  int
 *
 *  @ingroup attributes
 */
class VectorIntAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorIntAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorIntAttribute(std::vector<int> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        int  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<int> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<int>& i) {
        m_val = i;
    }

private:
    std::vector<int> m_val; ///< Attribute value
};

/**
 *  @class HepMC3::VectorIntAttribute
 *  @brief Attribute that holds a vector of integers of type  int
 *
 *  @ingroup attributes
 */
class VectorLongIntAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorLongIntAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorLongIntAttribute(std::vector<long int> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        long int  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<long int> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<long int>& i) {
        m_val = i;
    }

private:
    std::vector<long int> m_val; ///< Attribute value
};

/**
 *  @class HepMC3::VectorIntAttribute
 *  @brief Attribute that holds a vector of FPs of type  double
 *
 *  @ingroup attributes
 */
class VectorDoubleAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorDoubleAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorDoubleAttribute(std::vector<double> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att) override {
        double  datafoo;
        m_val.clear();
        std::stringstream datastream(att);
        while (datastream >> datafoo) m_val.push_back(datafoo);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=std::to_string(a);}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<double> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<double>& i) {
        m_val = i;
    }

private:
    std::vector<double> m_val; ///< Attribute value
};


/**
 *  @class HepMC3::VectorIntAttribute
 *  @brief Attribute that holds a vector of FPs of type  string
 *
 *  @ingroup attributes
 */
class VectorStringAttribute : public Attribute {
public:

    /** @brief Default constructor */
    VectorStringAttribute():Attribute(),m_val() {}

    /** @brief Constructor initializing attribute value */
    VectorStringAttribute(std::vector<std::string> val):Attribute(),m_val(val) {}

    /** @brief Implementation of Attribute::from_string */
    bool from_string(const string &att) override {
        size_t posb = att.find_first_not_of(' ');
        do {
           size_t pose = att.find_first_of(' ', posb);
           m_val.push_back(att.substr(posb, pose - posb));
           posb = att.find_first_not_of(' ', pose);
        } while (posb != std::string::npos);
        return true;
    }

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const  override{
        att.clear();
        for (auto a:  m_val) {if (att.length()) att+=" ";  att+=a;}
        return true;
    }

    /** @brief get the value associated to this Attribute. */
    std::vector<std::string> value() const {
        return m_val;
    }

    /** @brief set the value associated to this Attribute. */
    void set_value(const std::vector<std::string>& i) {
        m_val = i;
    }

private:
    std::vector<std::string> m_val; ///< Attribute value
};


} // namespace HepMC3

#endif
