#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/wx.h>
#include <wx/process.h>
#include <wx/statbmp.h>
#include <wx/image.h>
#include <wx/mimetype.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
#include "boost/regex.hpp"  
#include "boost/algorithm/string/regex.hpp"
namespace alg = boost::algorithm;
#include <fstream>
#include <string>
#include <vector>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/ClassMinerOptions.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Discretize.h"
#include "carmen/src/Dataset.h"
#include "carmen/src/Parser.h"
#include "PanelPlotter.h"
#include "Investigation.h"
using namespace boost;
using namespace std;
#include "DlgProgress.h"
#include "DlgPlot.h"

IMPLEMENT_CLASS( PlotDialog, wxDialog )

PlotDialog::PlotDialog( Investigation* investigation, Dataset* dataset){
	this->investigation = investigation;
    for(int i=0; i<(int)this->investigation->ga->identifiers.size(); i++){
		std::string identifier = this->investigation->ga->identifiers.at(i);
		std::string gene = this->investigation->ga->prop_for_identifier(identifier, this->investigation->gene_name_column);
		boost::algorithm::to_lower(gene);
		this->g_p.push_back( new std::pair<std::string, std::string>(gene, identifier) );
	}
    this->dataset = dataset;
	wxDialog::Create( NULL, wxID_ANY, wxString::FromAscii("Plot Expression Data"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	wxStaticText* lbl_probes = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Probes:") );
	this->txt_probes = new wxTextCtrl( this, ID_PLOT_TXT_PROBES, wxString::FromAscii(""), wxDefaultPosition, wxSize(150, 150), wxTE_MULTILINE );
    wxStaticText* lbl_gene = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Add Gene:") );
	
    this->txt_gene = new wxTextCtrl( this, ID_PLOT_TXT_GENE, wxString::FromAscii(""), wxDefaultPosition, wxSize(100, 20));
    this->btn_add_gene = new wxButton( this, ID_PLOT_BTN_ADD, wxString::FromAscii("&Add"), wxDefaultPosition, wxSize(50,20));

    wxStaticText* lbl_group_by = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Group By:") );
    wxArrayString ars_group_by, ars_limit_to, ars_is_or_is_not, ars_label_attribute;
    
    ars_group_by.Add(wxT("None"));
    ars_limit_to.Add(wxT("All"));
    ars_is_or_is_not.Add(wxT("is"));
    ars_is_or_is_not.Add(wxT("is not"));
    for( int i=0;i<(int)this->investigation->sa->attrib.size(); i++){
        ars_group_by.Add(wxString::FromAscii( this->investigation->sa->attrib.at(i).c_str() ));
        ars_limit_to.Add(wxString::FromAscii( this->investigation->sa->attrib.at(i).c_str() ));
        ars_label_attribute.Add(wxString::FromAscii( this->investigation->sa->attrib.at(i).c_str() ));
        limits.push_back(this->investigation->sa->attrib.at(i).c_str());
    }
    
	this->cho_group_by = new wxChoice( this, ID_PLOT_CHO_GROUP_BY, wxDefaultPosition, wxDefaultSize, ars_group_by );
	this->cho_group_by->SetSelection(0);

	wxStaticText* lbl_label_attribute = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Labels:") );
	this->cho_label_attribute = new wxChoice( this, ID_PLOT_CHO_LABEL_ATTRIBUTE, wxDefaultPosition, wxDefaultSize, ars_label_attribute );
	this->cho_label_attribute->SetSelection(0);

	wxStaticText* lbl_limit_display_to = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Limit to:") );
	this->cho_limit_to = new wxChoice( this, ID_PLOT_CHO_LIMIT_BY, wxDefaultPosition, wxDefaultSize, ars_limit_to );
	this->cho_limit_to->SetSelection(0);
    this->cho_is_or_is_not = new wxChoice( this, ID_PLOT_CHO_IS_OR_IS_NOT, wxDefaultPosition, wxSize(50, 20), ars_is_or_is_not);
    this->cho_is_or_is_not ->SetSelection(0);
    this->cho_limit_to_value = new wxChoice( this, ID_PLOT_CHO_LIMIT_TO_VALUE, wxDefaultPosition, wxSize(50, 20), ars_limit_to_value);
	this->cho_limit_to_value->SetSelection(0);
    this->cho_limit_to_value->Enable(false);
    this->cho_is_or_is_not->Enable(false);

    wxStaticText* lbl_sort = new wxStaticText( this, wxID_ANY, wxT("Sorted") );
	this->chk_sort = new wxCheckBox(this, ID_PLOT_CHK_SORT, wxString::FromAscii(""));

	wxStaticText* lbl_lines = new wxStaticText( this, wxID_ANY, wxT("Connect dots") );
	this->chk_lines = new wxCheckBox(this, ID_PLOT_CHK_LINES, wxString::FromAscii(""));

	wxStaticText* lbl_key_location = new wxStaticText( this, wxID_ANY, wxT("Key:") );
	this->ars_key_locations.Add(wxT("Top Right"));
	this->ars_key_locations.Add(wxT("Top Left"));
	this->ars_key_locations.Add(wxT("Lower Right"));
	this->ars_key_locations.Add(wxT("Lower Left"));
	this->cho_key_location = new wxChoice( this, ID_PLOT_KEY_LOCATION, wxDefaultPosition, wxSize(120,20), ars_key_locations);
	this->cho_key_location->SetSelection(0);

	this->ars_height.Add(wxT("600"));
	this->ars_height.Add(wxT("800"));
	this->ars_height.Add(wxT("1000"));
	this->ars_width.Add(wxT("800"));
	this->ars_width.Add(wxT("1000"));
	this->ars_width.Add(wxT("1200"));

	wxStaticText* lbl_height = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Height:") );
	this->cho_height = new wxChoice( this, ID_PLOT_CHO_HEIGHT, wxDefaultPosition, wxDefaultSize, ars_height );
	this->cho_height->SetSelection(0);

	wxStaticText* lbl_width = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Width:") );
	this->cho_width = new wxChoice( this, ID_PLOT_CHO_WIDTH, wxDefaultPosition, wxDefaultSize, ars_width );
	this->cho_width->SetSelection(0);

    this->btn_export = new wxButton( this, ID_PLOT_BTN_EXPORT, wxT("&Export"));
    this->btn_plot = new wxButton( this, ID_PLOT_BTN_PLOT, wxT("&Plot"));
    this->btn_close = new wxButton( this, wxID_CANCEL, wxT("&Close"));
    this->plotter = new Plotter(this, ID_PLOT_BMP, IMG_WIDTH, IMG_HEIGHT);

    Connect(ID_PLOT_BTN_ADD, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PlotDialog::OnClickAdd));
    Connect(ID_PLOT_BTN_PLOT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PlotDialog::OnClickPlot));
    Connect(ID_PLOT_BTN_EXPORT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PlotDialog::OnClickExport));
    Connect(wxID_CANCEL, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PlotDialog::OnClickClose));
    Connect(ID_PLOT_CHO_LIMIT_BY, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(PlotDialog::OnChangeAttribute));
    Connect(wxEVT_CHAR, wxKeyEventHandler(PlotDialog::OnChar));
    Connect(ID_PLOT_TXT_PROBES, wxEVT_CHAR, wxKeyEventHandler(PlotDialog::OnChar));
    Connect(ID_PLOT_TXT_GENE, wxEVT_CHAR, wxKeyEventHandler(PlotDialog::OnChar));

    wxBoxSizer* sizer_top = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* sizer_left = new wxBoxSizer(wxVERTICAL);
    sizer_top->Add(sizer_left, 0, wxALL | wxALIGN_LEFT | wxALIGN_TOP, BORDER_PXL);

    wxFlexGridSizer* sizer_options= new wxFlexGridSizer(10, 2, 0, 0);
    sizer_options->Add(lbl_probes, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(txt_probes, 0, FLAGS, BORDER_PXL);
    
    // add gene
    sizer_options->Add(lbl_gene, 0, FLAGS, BORDER_PXL);
    wxBoxSizer* sizer_gene = new wxBoxSizer(wxHORIZONTAL);
    sizer_gene->Add(txt_gene, 0, FLAGS, BORDER_PXL);
    sizer_gene->Add(btn_add_gene, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(sizer_gene, 0, FLAGS, BORDER_PXL);

    // group by
    sizer_options->Add(lbl_group_by, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(this->cho_group_by, 0, FLAGS, BORDER_PXL);
    
    // label by
    sizer_options->Add(lbl_label_attribute, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(this->cho_label_attribute, 0, FLAGS, BORDER_PXL);

    // limit display
    sizer_options->Add(lbl_limit_display_to, 0, FLAGS, BORDER_PXL);
    wxBoxSizer* sizer_limit_to = new wxBoxSizer(wxHORIZONTAL);
    sizer_limit_to->Add(this->cho_limit_to, 0, FLAGS, BORDER_NONE);
    sizer_limit_to->Add(this->cho_is_or_is_not, 0, FLAGS, BORDER_NONE);
    sizer_limit_to->Add(this->cho_limit_to_value, 0, FLAGS, BORDER_NONE);
    sizer_options->Add(sizer_limit_to, 0, FLAGS, BORDER_PXL);

    // sort and lines
    sizer_options->Add(lbl_sort, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(chk_sort, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(lbl_lines, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(chk_lines, 0, FLAGS, BORDER_PXL);

    // key location
    sizer_options->Add(lbl_key_location, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(this->cho_key_location, 0, FLAGS, BORDER_PXL);

    sizer_options->Add(lbl_height, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(this->cho_height, 0, FLAGS, BORDER_PXL);

    sizer_options->Add(lbl_width, 0, FLAGS, BORDER_PXL);
    sizer_options->Add(this->cho_width, 0, FLAGS, BORDER_PXL);

    sizer_left->Add(sizer_options, 0, FLAGS, BORDER_NONE);

    wxBoxSizer* sizer_btn = new wxBoxSizer(wxHORIZONTAL);
    sizer_btn->Add(btn_export, 0, FLAGS, BORDER_PXL);
    sizer_btn->Add(btn_plot, 0, FLAGS, BORDER_PXL);
    sizer_left->Add(sizer_btn, 0, FLAGS_RIGHT, BORDER_PXL);
    sizer_left->Add( new wxStaticText( this, wxID_ANY, wxString::FromAscii(""), wxDefaultPosition, wxSize(5,25) ));

    sizer_left->Add(btn_close, 0, FLAGS, BORDER_PXL);
    sizer_top->Add(this->plotter, 1, wxEXPAND);

    if(this->investigation->is_mac){
		// mac uses giant fonts by default; PC is small enough by default.
    	wxFont font(lbl_probes->GetFont());
    	font.SetPointSize(wxSMALL_FONT->GetPointSize());
		lbl_probes->SetFont(font);
		lbl_limit_display_to->SetFont(font);
		lbl_group_by->SetFont( font );
		lbl_sort->SetFont( font );
		lbl_lines->SetFont( font );
		lbl_height->SetFont( font );
		lbl_width->SetFont( font );
		this->cho_group_by->SetFont( font );
		this->btn_add_gene->SetFont( font );
		lbl_gene->SetFont( font );
		this->cho_limit_to->SetFont( font );
		this->chk_sort->SetFont( font );
		this->chk_lines->SetFont( font );
		this->btn_export->SetFont( font );
		this->btn_plot->SetFont( font );
		this->btn_close->SetFont( font );
		this->cho_key_location->SetFont(font);
		lbl_key_location->SetFont(font);
		this->cho_label_attribute->SetFont(font);
		lbl_label_attribute->SetFont(font);
	}

    this->SetSizer(sizer_top);
    this->SetAutoLayout(true);

    plotter_width=800;
    plotter_height=600;

	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	this->Refresh();
	this->Update();
	Centre();
	this->GetSize(&total_width, &total_height);
}


void PlotDialog::OnChar(wxKeyEvent& event){
	std::cout << event.GetKeyCode() << "\n";
	if (event.GetKeyCode() == WXK_ESCAPE){
		this->EndModal(wxCANCEL);
	}
	event.Skip();
};


void PlotDialog::OnClickPlot(wxCommandEvent& WXUNUSED(event)){
    wxString items = this->txt_probes->GetValue();
    if( (int)items.size()==0 ){
		wxMessageBox(wxString::FromAscii("ERROR: Please enter at least one probe."), _T("ERROR: no probes entered"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
    std::vector<std::string> probes, probes_ungroomed, cmd;
    std::string probes_string( items.ToAscii() );
    alg::split(probes_ungroomed, probes_string, alg::is_any_of("\n") );
    std::stringstream unknown_probes;
    for(int i=0; i<(int)probes_ungroomed.size(); i++){
        if( (int)probes_ungroomed.at(i).size()>0 ){
            if( this->investigation->ga->identifier2idx.find(probes_ungroomed.at(i)) == this->investigation->ga->identifier2idx.end() ){
                unknown_probes << probes_ungroomed.at(i) << " ";
            }
            else{
            	probes.push_back(probes_ungroomed.at(i));
            }
        }
    }
    if(unknown_probes.str().size()>0){
    	wxMessageBox(wxString::FromAscii(unknown_probes.str().c_str()), _T("ERROR: unknown probe"), wxOK | wxICON_EXCLAMATION, this);
    	if( probes.size()==0 )
    		return;
    }
    std::vector<std::string> labels, genes;
    std::string sample_limit = std::string("IDENTIFIER!NULL");
    if( this->cho_limit_to->GetSelection() != 0){
        std::stringstream ss;
        int sel = this->cho_limit_to->GetSelection()-1;
        ss << this->limits.at( sel );
        if( this->cho_is_or_is_not->GetSelection() == 0 )
            ss << "=";
        else
            ss << "!";
        ss << this->values.at( this->cho_limit_to_value->GetSelection() ) ;
        sample_limit = ss.str();
    }
    bool sort_by_value = false;
    if( this->chk_sort->IsChecked())
    	sort_by_value = true;
    std::string group_by("");
    int gb_idx = this->cho_group_by->GetCurrentSelection();
    if( gb_idx>0 )
    	group_by = this->investigation->sa->attrib[ gb_idx-1 ];
    int keyloc = this->cho_key_location->GetCurrentSelection();
    std::string label_attribute( this->investigation->sa->attrib[ this->cho_label_attribute->GetCurrentSelection() ] );

    try{
    	std::vector<std::vector<double>*> plot_values;
    	this->dataset->extract(plot_values, labels, genes, probes, sample_limit, group_by, sort_by_value, label_attribute );
    	for(int i=0; i<(int)genes.size(); i++){
    		std::stringstream ss;
    		ss << genes.at(i) << " " << probes.at(i);
    		genes.at(i) = ss.str();
    	}
    	this->plotter->set_key_location(keyloc);
    	this->plotter->set_plot_values(labels, genes, plot_values);
    	int height,width;
    	if( this->cho_height->GetSelection()==0 )
    		height=600;
    	else if( this->cho_height->GetSelection()==1 )
    		height=800;
    	else
    		height=1000;
    	if( this->cho_width->GetSelection()==0 )
    		width=800;
    	else if( this->cho_width->GetSelection()==1 )
    		width=1000;
    	else
    		width=1200;
    	this->plotter->DrawLines( this->chk_lines->IsChecked() );
    	this->plotter->SetDimensions(height,width);
    	this->plotter->Repaint();
    	for(int i=0; i<(int)plot_values.size(); i++){
    		if(plot_values.at(i) != NULL){
    			delete plot_values.at(i);
    		}
    	}
    	if( width != plotter_width || height != plotter_height){
    		total_width = total_width + ( width-plotter_width);
    		total_height = total_height + ( height-plotter_height);
    		plotter_width = width;
    		plotter_height = height;
    		this->SetSize(total_width, total_height);
    	}

    	this->Refresh();
    	this->Update();

    }
    catch(std::string err){
    	wxMessageBox(wxString::FromAscii(err.c_str()));
    }
}



void PlotDialog::ShowTextFile(std::string fn){
	if(this->investigation->is_mac){
		std::stringstream ss;
		ss << "open -t " << fn; // -t flag opens with default text editor
		wxExecute(wxString::FromAscii(ss.str().c_str()));
	}
	else{
		std::string loc_text = this->investigation->cp->get(std::string("External"), std::string("text_editor"));;
		std::stringstream ss;
		ss << loc_text << " " << fn;
		wxExecute(wxString::FromAscii(ss.str().c_str()));
	}
}


// recover code that requests file name

void PlotDialog::OnClickExport(wxCommandEvent& WXUNUSED(event)){
	wxString default_location = wxString::FromAscii(this->investigation->dir_results.c_str());
	wxString default_name = wxString::FromAscii("export.txt");
	wxString default_extension = wxString::FromAscii("txt");

	wxString fileloc = wxFileSelector(wxString::FromAscii("Save to:"), default_location, default_name, default_extension, wxString::FromAscii("*.*"), wxFD_SAVE, this);
	if( fileloc.size()==0 ){
		return;
	}
	std::string fn_out(fileloc.ToAscii());
	wxString items = this->txt_probes->GetValue();
    if( (int)items.size()==0 ){
		wxMessageBox(wxString::FromAscii("ERROR: Please enter at least one probe."), _T("ERROR: no probes entered"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
    std::vector<std::string> probes, cmd, probes_validated;
    std::string probes_raw( items.ToAscii() );
    alg::split(probes, probes_raw, alg::is_any_of("\n") );
    for(int i=0; i<(int)probes.size(); i++){
        if( (int)probes.at(i).size()>0 ){
            if( this->investigation->ga->identifier2idx.find(probes.at(i)) == this->investigation->ga->identifier2idx.end() ){
                std::string s("ERROR: unknown probe ");
                s += probes.at(i);
                wxMessageBox(wxString::FromAscii(s.c_str()), _T("ERROR: unknown probe"), wxOK | wxICON_EXCLAMATION, this);
		        return;
            }
            else{
            	probes_validated.push_back(probes.at(i));
            }
        }
    }

    std::vector<std::string> genes;
    std::string sample_limit = std::string("IDENTIFIER!NULL");
    if( this->cho_limit_to->GetSelection() != 0){
        std::stringstream ss;
        int sel = this->cho_limit_to->GetSelection()-1;
        ss << this->limits.at( sel );
        if( this->cho_is_or_is_not->GetSelection() == 0 )
            ss << "=";
        else
            ss << "!";
        ss << this->values.at( this->cho_limit_to_value->GetSelection() ) ;
        sample_limit = ss.str();
    }
    bool sort_by_value = false;
    if( this->chk_sort->IsChecked())
    	sort_by_value = true;
    std::string group_by("");
    int gb_idx = this->cho_group_by->GetCurrentSelection();
    if( gb_idx>0 )
    	group_by = this->investigation->sa->attrib[ gb_idx-1 ];

    try{
    	std::vector<std::vector<double>*> plot_values;
    	std::vector<std::string> labels;
    	this->dataset->extract(plot_values, labels, genes, probes_validated, sample_limit, group_by, sort_by_value, std::string("IDENTIFIER") );
    	this->investigation->write_dataset(fn_out, probes_validated, labels, 0);
    	this->ShowTextFile(fn_out);
    }
    catch(std::string err){
    	wxMessageBox(wxString::FromAscii(err.c_str()));
    }
}

void PlotDialog::OnChangeAttribute(wxCommandEvent& WXUNUSED(event)){
	this->cho_limit_to_value->Clear();
	this->values.clear();
    int idx_selection = this->cho_limit_to->GetSelection();
    if( idx_selection>0 ){
        this->cho_limit_to_value->Enable(true);
        this->cho_is_or_is_not->Enable(true);
        std::string attribute = this->investigation->sa->attrib.at( idx_selection-1 );
	    HASH_S_I seen;
	    std::string value;
    	
	    for(int i=0; i<(int)this->investigation->sa->identifiers.size(); i++){
		    value = this->investigation->sa->prop_for_identifier(this->investigation->sa->identifiers.at(i),attribute);
		    if(seen.find( value ) == seen.end()){
			    seen[value] = 1;
			    this->values.push_back(value);
		    }
	    }
	    std::sort(this->values.begin(), this->values.end());
	    for(int i=0; i<(int)this->values.size(); i++)
		    this->cho_limit_to_value->Append(wxString::FromAscii(this->values.at(i).c_str() ));
	    this->cho_limit_to_value->SetSelection(0);
    }
    else{
        this->cho_limit_to_value->Enable(false);
        this->cho_is_or_is_not->Enable(false);
    }
}

void PlotDialog::OnClickAdd(wxCommandEvent& WXUNUSED(event)){

    wxString search_w = this->txt_gene->GetValue();
    std::string searchterm(search_w.ToAscii());
	boost::algorithm::to_lower(searchterm);
    wxArrayString aProbes;
	for(int i=0; i<(int)this->g_p.size(); i++)
		if( this->g_p.at(i)->first.compare(searchterm)==0 )
            aProbes.Add( wxString::FromAscii( this->g_p.at(i)->second.c_str()) );
    if( aProbes.size()==0 ){
        wxMessageBox(wxString::FromAscii("ERROR: Gene symbol not found in this investigation."), _T("ERROR: gene not found"), wxOK | wxICON_EXCLAMATION, this);
		return;
    }
    wxString probe_id;
    if( aProbes.size()==1 ){
        probe_id = aProbes.Item(0);
    }
    else{
    	wxSingleChoiceDialog dialog(this, wxString::FromAscii("Choose a probe"), wxString::FromAscii("Choose a probe"), aProbes);
    	dialog.SetSelection(0);
    	if( dialog.ShowModal() == wxID_OK ){
    		probe_id = dialog.GetStringSelection();
    	}
        //probe_id = wxGetSingleChoice(wxString::FromAscii("Choose a probe"), wxString::FromAscii("Choose a probe"), aProbes);
    }
    if( probe_id.size()>0 ){
        wxString items = this->txt_probes->GetValue();
        if( items.size()==0 )
            this->txt_probes->SetValue( probe_id );
        else{
            std::string probes_raw( items.ToAscii() );
            std::vector<std::string> probes;
            alg::split(probes, probes_raw, alg::is_any_of("\n") );
            probes.push_back( std::string( probe_id.ToAscii() ) );
            std::string probe_list;
            stringhasher::join( probe_list, &probes, "\n" );
            this->txt_probes->SetValue( wxString::FromAscii(probe_list.c_str()) );
        }
    }
}

void PlotDialog::OnClickClose(wxCommandEvent& WXUNUSED(event)){
	for(int i=0; i<int(this->g_p.size()); i++){
		delete this->g_p.at(i);
	}
	this->AcceptAndClose();
}
