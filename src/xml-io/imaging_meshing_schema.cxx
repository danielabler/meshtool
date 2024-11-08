// Copyright (c) 2005-2014 Code Synthesis Tools CC
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis Tools CC gives
// permission to link this program with the Xerces-C++ library (or with
// modified versions of Xerces-C++ that use the same license as Xerces-C++),
// and distribute linked combinations including the two. You must obey
// the GNU General Public License version 2 in all respects for all of
// the code used other than Xerces-C++. If you modify this copy of the
// program, you may extend this exception to your version of the program,
// but you are not obligated to do so. If you do not wish to do so, delete
// this exception statement from your version.
//
// Furthermore, Code Synthesis Tools CC makes a special exception for
// the Free/Libre and Open Source Software (FLOSS) which is described
// in the accompanying FLOSSE file.
//

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/pre.hxx>

#include "imaging_meshing_schema.hxx"

// MeshTool_Config_t
// 

const MeshTool_Config_t::Operation_Image2Mesh_optional& MeshTool_Config_t::
Operation_Image2Mesh () const
{
  return this->Operation_Image2Mesh_;
}

MeshTool_Config_t::Operation_Image2Mesh_optional& MeshTool_Config_t::
Operation_Image2Mesh ()
{
  return this->Operation_Image2Mesh_;
}

void MeshTool_Config_t::
Operation_Image2Mesh (const Operation_Image2Mesh_type& x)
{
  this->Operation_Image2Mesh_.set (x);
}

void MeshTool_Config_t::
Operation_Image2Mesh (const Operation_Image2Mesh_optional& x)
{
  this->Operation_Image2Mesh_ = x;
}

void MeshTool_Config_t::
Operation_Image2Mesh (::std::auto_ptr< Operation_Image2Mesh_type > x)
{
  this->Operation_Image2Mesh_.set (x);
}

const MeshTool_Config_t::Operation_Remesh_optional& MeshTool_Config_t::
Operation_Remesh () const
{
  return this->Operation_Remesh_;
}

MeshTool_Config_t::Operation_Remesh_optional& MeshTool_Config_t::
Operation_Remesh ()
{
  return this->Operation_Remesh_;
}

void MeshTool_Config_t::
Operation_Remesh (const Operation_Remesh_type& x)
{
  this->Operation_Remesh_.set (x);
}

void MeshTool_Config_t::
Operation_Remesh (const Operation_Remesh_optional& x)
{
  this->Operation_Remesh_ = x;
}

void MeshTool_Config_t::
Operation_Remesh (::std::auto_ptr< Operation_Remesh_type > x)
{
  this->Operation_Remesh_.set (x);
}


// Operation_Image2Mesh_t
// 

const Operation_Image2Mesh_t::MeshCriteria_global_type& Operation_Image2Mesh_t::
MeshCriteria_global () const
{
  return this->MeshCriteria_global_.get ();
}

Operation_Image2Mesh_t::MeshCriteria_global_type& Operation_Image2Mesh_t::
MeshCriteria_global ()
{
  return this->MeshCriteria_global_.get ();
}

void Operation_Image2Mesh_t::
MeshCriteria_global (const MeshCriteria_global_type& x)
{
  this->MeshCriteria_global_.set (x);
}

void Operation_Image2Mesh_t::
MeshCriteria_global (::std::auto_ptr< MeshCriteria_global_type > x)
{
  this->MeshCriteria_global_.set (x);
}

const Operation_Image2Mesh_t::MeshCriteria_SubDomain_sequence& Operation_Image2Mesh_t::
MeshCriteria_SubDomain () const
{
  return this->MeshCriteria_SubDomain_;
}

Operation_Image2Mesh_t::MeshCriteria_SubDomain_sequence& Operation_Image2Mesh_t::
MeshCriteria_SubDomain ()
{
  return this->MeshCriteria_SubDomain_;
}

void Operation_Image2Mesh_t::
MeshCriteria_SubDomain (const MeshCriteria_SubDomain_sequence& s)
{
  this->MeshCriteria_SubDomain_ = s;
}

const Operation_Image2Mesh_t::path_to_input_file_type& Operation_Image2Mesh_t::
path_to_input_file () const
{
  return this->path_to_input_file_.get ();
}

Operation_Image2Mesh_t::path_to_input_file_type& Operation_Image2Mesh_t::
path_to_input_file ()
{
  return this->path_to_input_file_.get ();
}

void Operation_Image2Mesh_t::
path_to_input_file (const path_to_input_file_type& x)
{
  this->path_to_input_file_.set (x);
}

void Operation_Image2Mesh_t::
path_to_input_file (::std::auto_ptr< path_to_input_file_type > x)
{
  this->path_to_input_file_.set (x);
}

const Operation_Image2Mesh_t::path_to_output_file_type& Operation_Image2Mesh_t::
path_to_output_file () const
{
  return this->path_to_output_file_.get ();
}

Operation_Image2Mesh_t::path_to_output_file_type& Operation_Image2Mesh_t::
path_to_output_file ()
{
  return this->path_to_output_file_.get ();
}

void Operation_Image2Mesh_t::
path_to_output_file (const path_to_output_file_type& x)
{
  this->path_to_output_file_.set (x);
}

void Operation_Image2Mesh_t::
path_to_output_file (::std::auto_ptr< path_to_output_file_type > x)
{
  this->path_to_output_file_.set (x);
}


// MeshCriteria_global_t
// 

const MeshCriteria_global_t::facet_angle_type& MeshCriteria_global_t::
facet_angle () const
{
  return this->facet_angle_.get ();
}

MeshCriteria_global_t::facet_angle_type& MeshCriteria_global_t::
facet_angle ()
{
  return this->facet_angle_.get ();
}

void MeshCriteria_global_t::
facet_angle (const facet_angle_type& x)
{
  this->facet_angle_.set (x);
}

const MeshCriteria_global_t::facet_size_type& MeshCriteria_global_t::
facet_size () const
{
  return this->facet_size_.get ();
}

MeshCriteria_global_t::facet_size_type& MeshCriteria_global_t::
facet_size ()
{
  return this->facet_size_.get ();
}

void MeshCriteria_global_t::
facet_size (const facet_size_type& x)
{
  this->facet_size_.set (x);
}

const MeshCriteria_global_t::facet_distance_type& MeshCriteria_global_t::
facet_distance () const
{
  return this->facet_distance_.get ();
}

MeshCriteria_global_t::facet_distance_type& MeshCriteria_global_t::
facet_distance ()
{
  return this->facet_distance_.get ();
}

void MeshCriteria_global_t::
facet_distance (const facet_distance_type& x)
{
  this->facet_distance_.set (x);
}

const MeshCriteria_global_t::cell_radius_edge_ratio_type& MeshCriteria_global_t::
cell_radius_edge_ratio () const
{
  return this->cell_radius_edge_ratio_.get ();
}

MeshCriteria_global_t::cell_radius_edge_ratio_type& MeshCriteria_global_t::
cell_radius_edge_ratio ()
{
  return this->cell_radius_edge_ratio_.get ();
}

void MeshCriteria_global_t::
cell_radius_edge_ratio (const cell_radius_edge_ratio_type& x)
{
  this->cell_radius_edge_ratio_.set (x);
}

const MeshCriteria_global_t::cell_size_type& MeshCriteria_global_t::
cell_size () const
{
  return this->cell_size_.get ();
}

MeshCriteria_global_t::cell_size_type& MeshCriteria_global_t::
cell_size ()
{
  return this->cell_size_.get ();
}

void MeshCriteria_global_t::
cell_size (const cell_size_type& x)
{
  this->cell_size_.set (x);
}


// MeshCriteria_SubDomain_t
// 

const MeshCriteria_SubDomain_t::domain_id_type& MeshCriteria_SubDomain_t::
domain_id () const
{
  return this->domain_id_.get ();
}

MeshCriteria_SubDomain_t::domain_id_type& MeshCriteria_SubDomain_t::
domain_id ()
{
  return this->domain_id_.get ();
}

void MeshCriteria_SubDomain_t::
domain_id (const domain_id_type& x)
{
  this->domain_id_.set (x);
}

const MeshCriteria_SubDomain_t::cell_size_type& MeshCriteria_SubDomain_t::
cell_size () const
{
  return this->cell_size_.get ();
}

MeshCriteria_SubDomain_t::cell_size_type& MeshCriteria_SubDomain_t::
cell_size ()
{
  return this->cell_size_.get ();
}

void MeshCriteria_SubDomain_t::
cell_size (const cell_size_type& x)
{
  this->cell_size_.set (x);
}

const MeshCriteria_SubDomain_t::dimension_type& MeshCriteria_SubDomain_t::
dimension () const
{
  return this->dimension_.get ();
}

MeshCriteria_SubDomain_t::dimension_type& MeshCriteria_SubDomain_t::
dimension ()
{
  return this->dimension_.get ();
}

void MeshCriteria_SubDomain_t::
dimension (const dimension_type& x)
{
  this->dimension_.set (x);
}

MeshCriteria_SubDomain_t::dimension_type MeshCriteria_SubDomain_t::
dimension_default_value ()
{
  return dimension_type (3);
}

const MeshCriteria_SubDomain_t::name_optional& MeshCriteria_SubDomain_t::
name () const
{
  return this->name_;
}

MeshCriteria_SubDomain_t::name_optional& MeshCriteria_SubDomain_t::
name ()
{
  return this->name_;
}

void MeshCriteria_SubDomain_t::
name (const name_type& x)
{
  this->name_.set (x);
}

void MeshCriteria_SubDomain_t::
name (const name_optional& x)
{
  this->name_ = x;
}

void MeshCriteria_SubDomain_t::
name (::std::auto_ptr< name_type > x)
{
  this->name_.set (x);
}


// Operation_Remesh_t
// 

const Operation_Remesh_t::spacing_xyz_type& Operation_Remesh_t::
spacing_xyz () const
{
  return this->spacing_xyz_.get ();
}

Operation_Remesh_t::spacing_xyz_type& Operation_Remesh_t::
spacing_xyz ()
{
  return this->spacing_xyz_.get ();
}

void Operation_Remesh_t::
spacing_xyz (const spacing_xyz_type& x)
{
  this->spacing_xyz_.set (x);
}


#include <xsd/cxx/xml/dom/parsing-source.hxx>

// MeshTool_Config_t
//

MeshTool_Config_t::
MeshTool_Config_t ()
: ::xml_schema::type (),
  Operation_Image2Mesh_ (this),
  Operation_Remesh_ (this)
{
}

MeshTool_Config_t::
MeshTool_Config_t (const MeshTool_Config_t& x,
                   ::xml_schema::flags f,
                   ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  Operation_Image2Mesh_ (x.Operation_Image2Mesh_, f, this),
  Operation_Remesh_ (x.Operation_Remesh_, f, this)
{
}

MeshTool_Config_t::
MeshTool_Config_t (const ::xercesc::DOMElement& e,
                   ::xml_schema::flags f,
                   ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  Operation_Image2Mesh_ (this),
  Operation_Remesh_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void MeshTool_Config_t::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // Operation_Image2Mesh
    //
    if (n.name () == "Operation_Image2Mesh" && n.namespace_ ().empty ())
    {
      ::std::auto_ptr< Operation_Image2Mesh_type > r (
        Operation_Image2Mesh_traits::create (i, f, this));

      if (!this->Operation_Image2Mesh_)
      {
        this->Operation_Image2Mesh_.set (r);
        continue;
      }
    }

    // Operation_Remesh
    //
    if (n.name () == "Operation_Remesh" && n.namespace_ ().empty ())
    {
      ::std::auto_ptr< Operation_Remesh_type > r (
        Operation_Remesh_traits::create (i, f, this));

      if (!this->Operation_Remesh_)
      {
        this->Operation_Remesh_.set (r);
        continue;
      }
    }

    break;
  }
}

MeshTool_Config_t* MeshTool_Config_t::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class MeshTool_Config_t (*this, f, c);
}

MeshTool_Config_t& MeshTool_Config_t::
operator= (const MeshTool_Config_t& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->Operation_Image2Mesh_ = x.Operation_Image2Mesh_;
    this->Operation_Remesh_ = x.Operation_Remesh_;
  }

  return *this;
}

MeshTool_Config_t::
~MeshTool_Config_t ()
{
}

// Operation_Image2Mesh_t
//

Operation_Image2Mesh_t::
Operation_Image2Mesh_t (const MeshCriteria_global_type& MeshCriteria_global,
                        const path_to_input_file_type& path_to_input_file,
                        const path_to_output_file_type& path_to_output_file)
: ::xml_schema::type (),
  MeshCriteria_global_ (MeshCriteria_global, this),
  MeshCriteria_SubDomain_ (this),
  path_to_input_file_ (path_to_input_file, this),
  path_to_output_file_ (path_to_output_file, this)
{
}

Operation_Image2Mesh_t::
Operation_Image2Mesh_t (::std::auto_ptr< MeshCriteria_global_type > MeshCriteria_global,
                        const path_to_input_file_type& path_to_input_file,
                        const path_to_output_file_type& path_to_output_file)
: ::xml_schema::type (),
  MeshCriteria_global_ (MeshCriteria_global, this),
  MeshCriteria_SubDomain_ (this),
  path_to_input_file_ (path_to_input_file, this),
  path_to_output_file_ (path_to_output_file, this)
{
}

Operation_Image2Mesh_t::
Operation_Image2Mesh_t (const Operation_Image2Mesh_t& x,
                        ::xml_schema::flags f,
                        ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  MeshCriteria_global_ (x.MeshCriteria_global_, f, this),
  MeshCriteria_SubDomain_ (x.MeshCriteria_SubDomain_, f, this),
  path_to_input_file_ (x.path_to_input_file_, f, this),
  path_to_output_file_ (x.path_to_output_file_, f, this)
{
}

Operation_Image2Mesh_t::
Operation_Image2Mesh_t (const ::xercesc::DOMElement& e,
                        ::xml_schema::flags f,
                        ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  MeshCriteria_global_ (this),
  MeshCriteria_SubDomain_ (this),
  path_to_input_file_ (this),
  path_to_output_file_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, true);
    this->parse (p, f);
  }
}

void Operation_Image2Mesh_t::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // MeshCriteria_global
    //
    if (n.name () == "MeshCriteria_global" && n.namespace_ ().empty ())
    {
      ::std::auto_ptr< MeshCriteria_global_type > r (
        MeshCriteria_global_traits::create (i, f, this));

      if (!MeshCriteria_global_.present ())
      {
        this->MeshCriteria_global_.set (r);
        continue;
      }
    }

    // MeshCriteria_SubDomain
    //
    if (n.name () == "MeshCriteria_SubDomain" && n.namespace_ ().empty ())
    {
      ::std::auto_ptr< MeshCriteria_SubDomain_type > r (
        MeshCriteria_SubDomain_traits::create (i, f, this));

      this->MeshCriteria_SubDomain_.push_back (r);
      continue;
    }

    break;
  }

  if (!MeshCriteria_global_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "MeshCriteria_global",
      "");
  }

  while (p.more_attributes ())
  {
    const ::xercesc::DOMAttr& i (p.next_attribute ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    if (n.name () == "path_to_input_file" && n.namespace_ ().empty ())
    {
      this->path_to_input_file_.set (path_to_input_file_traits::create (i, f, this));
      continue;
    }

    if (n.name () == "path_to_output_file" && n.namespace_ ().empty ())
    {
      this->path_to_output_file_.set (path_to_output_file_traits::create (i, f, this));
      continue;
    }
  }

  if (!path_to_input_file_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "path_to_input_file",
      "");
  }

  if (!path_to_output_file_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "path_to_output_file",
      "");
  }
}

Operation_Image2Mesh_t* Operation_Image2Mesh_t::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class Operation_Image2Mesh_t (*this, f, c);
}

Operation_Image2Mesh_t& Operation_Image2Mesh_t::
operator= (const Operation_Image2Mesh_t& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->MeshCriteria_global_ = x.MeshCriteria_global_;
    this->MeshCriteria_SubDomain_ = x.MeshCriteria_SubDomain_;
    this->path_to_input_file_ = x.path_to_input_file_;
    this->path_to_output_file_ = x.path_to_output_file_;
  }

  return *this;
}

Operation_Image2Mesh_t::
~Operation_Image2Mesh_t ()
{
}

// MeshCriteria_global_t
//

MeshCriteria_global_t::
MeshCriteria_global_t (const facet_angle_type& facet_angle,
                       const facet_size_type& facet_size,
                       const facet_distance_type& facet_distance,
                       const cell_radius_edge_ratio_type& cell_radius_edge_ratio,
                       const cell_size_type& cell_size)
: ::xml_schema::type (),
  facet_angle_ (facet_angle, this),
  facet_size_ (facet_size, this),
  facet_distance_ (facet_distance, this),
  cell_radius_edge_ratio_ (cell_radius_edge_ratio, this),
  cell_size_ (cell_size, this)
{
}

MeshCriteria_global_t::
MeshCriteria_global_t (const MeshCriteria_global_t& x,
                       ::xml_schema::flags f,
                       ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  facet_angle_ (x.facet_angle_, f, this),
  facet_size_ (x.facet_size_, f, this),
  facet_distance_ (x.facet_distance_, f, this),
  cell_radius_edge_ratio_ (x.cell_radius_edge_ratio_, f, this),
  cell_size_ (x.cell_size_, f, this)
{
}

MeshCriteria_global_t::
MeshCriteria_global_t (const ::xercesc::DOMElement& e,
                       ::xml_schema::flags f,
                       ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  facet_angle_ (this),
  facet_size_ (this),
  facet_distance_ (this),
  cell_radius_edge_ratio_ (this),
  cell_size_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, false, false, true);
    this->parse (p, f);
  }
}

void MeshCriteria_global_t::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  while (p.more_attributes ())
  {
    const ::xercesc::DOMAttr& i (p.next_attribute ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    if (n.name () == "facet_angle" && n.namespace_ ().empty ())
    {
      this->facet_angle_.set (facet_angle_traits::create (i, f, this));
      continue;
    }

    if (n.name () == "facet_size" && n.namespace_ ().empty ())
    {
      this->facet_size_.set (facet_size_traits::create (i, f, this));
      continue;
    }

    if (n.name () == "facet_distance" && n.namespace_ ().empty ())
    {
      this->facet_distance_.set (facet_distance_traits::create (i, f, this));
      continue;
    }

    if (n.name () == "cell_radius_edge_ratio" && n.namespace_ ().empty ())
    {
      this->cell_radius_edge_ratio_.set (cell_radius_edge_ratio_traits::create (i, f, this));
      continue;
    }

    if (n.name () == "cell_size" && n.namespace_ ().empty ())
    {
      this->cell_size_.set (cell_size_traits::create (i, f, this));
      continue;
    }
  }

  if (!facet_angle_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "facet_angle",
      "");
  }

  if (!facet_size_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "facet_size",
      "");
  }

  if (!facet_distance_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "facet_distance",
      "");
  }

  if (!cell_radius_edge_ratio_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "cell_radius_edge_ratio",
      "");
  }

  if (!cell_size_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "cell_size",
      "");
  }
}

MeshCriteria_global_t* MeshCriteria_global_t::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class MeshCriteria_global_t (*this, f, c);
}

MeshCriteria_global_t& MeshCriteria_global_t::
operator= (const MeshCriteria_global_t& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->facet_angle_ = x.facet_angle_;
    this->facet_size_ = x.facet_size_;
    this->facet_distance_ = x.facet_distance_;
    this->cell_radius_edge_ratio_ = x.cell_radius_edge_ratio_;
    this->cell_size_ = x.cell_size_;
  }

  return *this;
}

MeshCriteria_global_t::
~MeshCriteria_global_t ()
{
}

// MeshCriteria_SubDomain_t
//

MeshCriteria_SubDomain_t::
MeshCriteria_SubDomain_t (const domain_id_type& domain_id,
                          const cell_size_type& cell_size)
: ::xml_schema::type (),
  domain_id_ (domain_id, this),
  cell_size_ (cell_size, this),
  dimension_ (dimension_default_value (), this),
  name_ (this)
{
}

MeshCriteria_SubDomain_t::
MeshCriteria_SubDomain_t (const MeshCriteria_SubDomain_t& x,
                          ::xml_schema::flags f,
                          ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  domain_id_ (x.domain_id_, f, this),
  cell_size_ (x.cell_size_, f, this),
  dimension_ (x.dimension_, f, this),
  name_ (x.name_, f, this)
{
}

MeshCriteria_SubDomain_t::
MeshCriteria_SubDomain_t (const ::xercesc::DOMElement& e,
                          ::xml_schema::flags f,
                          ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  domain_id_ (this),
  cell_size_ (this),
  dimension_ (this),
  name_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, false, false, true);
    this->parse (p, f);
  }
}

void MeshCriteria_SubDomain_t::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  while (p.more_attributes ())
  {
    const ::xercesc::DOMAttr& i (p.next_attribute ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    if (n.name () == "domain_id" && n.namespace_ ().empty ())
    {
      this->domain_id_.set (domain_id_traits::create (i, f, this));
      continue;
    }

    if (n.name () == "cell_size" && n.namespace_ ().empty ())
    {
      this->cell_size_.set (cell_size_traits::create (i, f, this));
      continue;
    }

    if (n.name () == "dimension" && n.namespace_ ().empty ())
    {
      this->dimension_.set (dimension_traits::create (i, f, this));
      continue;
    }

    if (n.name () == "name" && n.namespace_ ().empty ())
    {
      this->name_.set (name_traits::create (i, f, this));
      continue;
    }
  }

  if (!domain_id_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "domain_id",
      "");
  }

  if (!cell_size_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "cell_size",
      "");
  }

  if (!dimension_.present ())
  {
    this->dimension_.set (dimension_default_value ());
  }
}

MeshCriteria_SubDomain_t* MeshCriteria_SubDomain_t::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class MeshCriteria_SubDomain_t (*this, f, c);
}

MeshCriteria_SubDomain_t& MeshCriteria_SubDomain_t::
operator= (const MeshCriteria_SubDomain_t& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->domain_id_ = x.domain_id_;
    this->cell_size_ = x.cell_size_;
    this->dimension_ = x.dimension_;
    this->name_ = x.name_;
  }

  return *this;
}

MeshCriteria_SubDomain_t::
~MeshCriteria_SubDomain_t ()
{
}

// Operation_Remesh_t
//

Operation_Remesh_t::
Operation_Remesh_t (const MeshCriteria_global_type& MeshCriteria_global,
                    const path_to_input_file_type& path_to_input_file,
                    const path_to_output_file_type& path_to_output_file,
                    const spacing_xyz_type& spacing_xyz)
: ::Operation_Image2Mesh_t (MeshCriteria_global,
                            path_to_input_file,
                            path_to_output_file),
  spacing_xyz_ (spacing_xyz, this)
{
}

Operation_Remesh_t::
Operation_Remesh_t (::std::auto_ptr< MeshCriteria_global_type > MeshCriteria_global,
                    const path_to_input_file_type& path_to_input_file,
                    const path_to_output_file_type& path_to_output_file,
                    const spacing_xyz_type& spacing_xyz)
: ::Operation_Image2Mesh_t (MeshCriteria_global,
                            path_to_input_file,
                            path_to_output_file),
  spacing_xyz_ (spacing_xyz, this)
{
}

Operation_Remesh_t::
Operation_Remesh_t (const Operation_Remesh_t& x,
                    ::xml_schema::flags f,
                    ::xml_schema::container* c)
: ::Operation_Image2Mesh_t (x, f, c),
  spacing_xyz_ (x.spacing_xyz_, f, this)
{
}

Operation_Remesh_t::
Operation_Remesh_t (const ::xercesc::DOMElement& e,
                    ::xml_schema::flags f,
                    ::xml_schema::container* c)
: ::Operation_Image2Mesh_t (e, f | ::xml_schema::flags::base, c),
  spacing_xyz_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, true);
    this->parse (p, f);
  }
}

void Operation_Remesh_t::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  this->::Operation_Image2Mesh_t::parse (p, f);

  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // spacing_xyz
    //
    if (n.name () == "spacing_xyz" && n.namespace_ ().empty ())
    {
      if (!spacing_xyz_.present ())
      {
        this->spacing_xyz_.set (spacing_xyz_traits::create (i, f, this));
        continue;
      }
    }

    break;
  }

  if (!spacing_xyz_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "spacing_xyz",
      "");
  }
}

Operation_Remesh_t* Operation_Remesh_t::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class Operation_Remesh_t (*this, f, c);
}

Operation_Remesh_t& Operation_Remesh_t::
operator= (const Operation_Remesh_t& x)
{
  if (this != &x)
  {
    static_cast< ::Operation_Image2Mesh_t& > (*this) = x;
    this->spacing_xyz_ = x.spacing_xyz_;
  }

  return *this;
}

Operation_Remesh_t::
~Operation_Remesh_t ()
{
}

#include <istream>
#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/tree/error-handler.hxx>

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (const ::std::string& u,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::tree::error_handler< char > h;

  ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

  return ::std::auto_ptr< ::MeshTool_Config_t > (
    ::MeshTool_Config (
      d, f | ::xml_schema::flags::own_dom, p));
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (const ::std::string& u,
                 ::xml_schema::error_handler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::auto_ptr< ::MeshTool_Config_t > (
    ::MeshTool_Config (
      d, f | ::xml_schema::flags::own_dom, p));
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (const ::std::string& u,
                 ::xercesc::DOMErrorHandler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::auto_ptr< ::MeshTool_Config_t > (
    ::MeshTool_Config (
      d, f | ::xml_schema::flags::own_dom, p));
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::std::istream& is,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::MeshTool_Config (isrc, f, p);
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::std::istream& is,
                 ::xml_schema::error_handler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::MeshTool_Config (isrc, h, f, p);
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::std::istream& is,
                 ::xercesc::DOMErrorHandler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::MeshTool_Config (isrc, h, f, p);
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::std::istream& is,
                 const ::std::string& sid,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::MeshTool_Config (isrc, f, p);
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::std::istream& is,
                 const ::std::string& sid,
                 ::xml_schema::error_handler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::MeshTool_Config (isrc, h, f, p);
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::std::istream& is,
                 const ::std::string& sid,
                 ::xercesc::DOMErrorHandler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::MeshTool_Config (isrc, h, f, p);
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::xercesc::InputSource& i,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xsd::cxx::tree::error_handler< char > h;

  ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

  return ::std::auto_ptr< ::MeshTool_Config_t > (
    ::MeshTool_Config (
      d, f | ::xml_schema::flags::own_dom, p));
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::xercesc::InputSource& i,
                 ::xml_schema::error_handler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::auto_ptr< ::MeshTool_Config_t > (
    ::MeshTool_Config (
      d, f | ::xml_schema::flags::own_dom, p));
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::xercesc::InputSource& i,
                 ::xercesc::DOMErrorHandler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::auto_ptr< ::MeshTool_Config_t > (
    ::MeshTool_Config (
      d, f | ::xml_schema::flags::own_dom, p));
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (const ::xercesc::DOMDocument& doc,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
{
  if (f & ::xml_schema::flags::keep_dom)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      static_cast< ::xercesc::DOMDocument* > (doc.cloneNode (true)));

    return ::std::auto_ptr< ::MeshTool_Config_t > (
      ::MeshTool_Config (
        d, f | ::xml_schema::flags::own_dom, p));
  }

  const ::xercesc::DOMElement& e (*doc.getDocumentElement ());
  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (n.name () == "MeshTool_Config" &&
      n.namespace_ () == "")
  {
    ::std::auto_ptr< ::MeshTool_Config_t > r (
      ::xsd::cxx::tree::traits< ::MeshTool_Config_t, char >::create (
        e, f, 0));
    return r;
  }

  throw ::xsd::cxx::tree::unexpected_element < char > (
    n.name (),
    n.namespace_ (),
    "MeshTool_Config",
    "");
}

::std::auto_ptr< ::MeshTool_Config_t >
MeshTool_Config (::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties&)
{
  ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > c (
    ((f & ::xml_schema::flags::keep_dom) &&
     !(f & ::xml_schema::flags::own_dom))
    ? static_cast< ::xercesc::DOMDocument* > (d->cloneNode (true))
    : 0);

  ::xercesc::DOMDocument& doc (c.get () ? *c : *d);
  const ::xercesc::DOMElement& e (*doc.getDocumentElement ());

  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (f & ::xml_schema::flags::keep_dom)
    doc.setUserData (::xml_schema::dom::tree_node_key,
                     (c.get () ? &c : &d),
                     0);

  if (n.name () == "MeshTool_Config" &&
      n.namespace_ () == "")
  {
    ::std::auto_ptr< ::MeshTool_Config_t > r (
      ::xsd::cxx::tree::traits< ::MeshTool_Config_t, char >::create (
        e, f, 0));
    return r;
  }

  throw ::xsd::cxx::tree::unexpected_element < char > (
    n.name (),
    n.namespace_ (),
    "MeshTool_Config",
    "");
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

