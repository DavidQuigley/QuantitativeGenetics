#define _CRT_SECURE_NO_DEPRECATE 1
#include <signal.h>
#include <string>
#include <vector>
#include <sstream>
#include <wx/wx.h>
#include <wx/timer.h>
#include <wx/dialog.h>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#include <wx/sizer.h>
#include <wx/filename.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/listbox.h>
#include <wx/process.h>
#include <wx/timer.h>
#include <wx/mimetype.h>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"

namespace alg = boost::algorithm;
#include "DlgProgress.h"

IMPLEMENT_CLASS( ProgressDialog, wxDialog )

long ProgressDialog::popen2(std::string cmd, int &fd){
	// This code is only called from the macintosh; the PC side uses wxExecute
	// This is because wxExecute seems only to know how to call .app packages
	// on the OS X side.
    
#ifdef WIN32
	return -1;  // This code would never be called from windows
#else
	int	writepipe[2] = {-1,-1},	/* parent -> child */
		readpipe [2] = {-1,-1};	/* child -> parent */
	long childpid;
	writepipe[0] = -1;
	if (pipe(readpipe) < 0 || pipe(writepipe) < 0)
	    return -1;
	#define	PARENT_READ readpipe[0]
	#define	CHILD_WRITE readpipe[1]
	#define CHILD_READ writepipe[0]
	#define PARENT_WRITE writepipe[1]

	childpid = fork();
	if (childpid < 0)
	    return childpid; // Fatal error, cannot fork
	else if (childpid == 0){
		// in child
	    close(PARENT_WRITE);
	    close(PARENT_READ);
	    dup2(CHILD_READ,  0);
	    close(CHILD_READ);
	    dup2(CHILD_WRITE, 1);
	    close(CHILD_WRITE);
        execl("/bin/sh", "sh", "-c", cmd.c_str(), NULL);
	    // should not get to this point, as execl only returns on fatal error
	    std::stringstream ss;
	    ss << "Fatal error in command: " << cmd.c_str();
	    wxMessageBox(wxString::FromAscii(ss.str().c_str()));
	    perror("execl");
	}
	else{
		// in parent
		close(CHILD_READ);
		close(CHILD_WRITE);
		fd = PARENT_READ;
		return childpid;
	}
	return -1;
#endif
}


ProgressDialog::ProgressDialog( wxWindow* parent, std::vector<std::string>& cmd, Investigation* investigation, bool show_progress ){
	
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( parent, wxID_ANY, wxT("Executing request..."));
	if(this->investigation->is_mac)
		this->SetWindowVariant(wxWINDOW_VARIANT_MINI);
	wxStaticText* lbl_progress = new wxStaticText(this, wxID_ANY, wxString::FromAscii("Progress:"));
    this->lbx_progress = new wxListBox(this, ID_LBX_PROGRESS, wxDefaultPosition, wxSize(450, 300) );
    if(this->investigation->is_mac){
        // mac uses giant fonts by default; PC is better.
    	wxFont font(this->lbx_progress->GetFont());
    	font.SetPointSize(wxSMALL_FONT->GetPointSize());
        this->lbx_progress->SetFont( font );
    }
	wxButton* btn_cancel = new wxButton( this, ID_BTN_CANCEL, wxT("&Cancel"));
	this->cancelled = false;
	this->error = "";

	wxBoxSizer* grid = new wxBoxSizer(wxVERTICAL);
    grid->Add( lbl_progress, 0, FLAGS, BORDER_PXL);
	grid->Add( this->lbx_progress, 0, FLAGS, BORDER_PXL);
	wxBoxSizer* sizer_btn = new wxBoxSizer(wxHORIZONTAL);
	sizer_btn->Add(btn_cancel, 0, FLAGS_RIGHT, BORDER_PXL );
	grid->Add( sizer_btn, 0, FLAGS, BORDER_PXL);
	this->SetSizer(grid);
	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	Centre(); 

	Connect(ID_BTN_CANCEL, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(ProgressDialog::OnClickCancel));
	Connect(ID_PROCESS, wxEVT_END_PROCESS, wxProcessEventHandler(ProgressDialog::OnProcessEnded));
	Connect(ID_PROCESS, wxEVT_IDLE, wxIdleEventHandler(ProgressDialog::OnIdle) );
	Connect(ID_TIMER, wxEVT_TIMER, wxTimerEventHandler(ProgressDialog::tick) );
	std::stringstream ss;

	if( cmd.at(0).compare("BATCH") == 0 ){
		if( this->investigation->is_mac )
			ss << "sh "; // Macintosh needs prefix of sh; PC will run Batch file without it.
	}
	else{
		// if Mac, call executable that is in the application bundle
		// if PC, call from homedir (Carmen directory) and append ".exe "
		ss << "\"";
		if( this->investigation->is_mac )
			ss << this->investigation->cp->get("Internal", "homedir") << "/Carmen.app/Contents/MacOS/";
		else
			ss << this->investigation->cp->get("Internal", "homedir") << "\\";
		ss << cmd.at(0);
		if( this->investigation->is_mac ){
            // check that the file we're about to run is executable
            std::stringstream ss_fn; // ss itself has a quote prepended to it, that hoses us.
            ss_fn << this->investigation->cp->get("Internal", "homedir") << "/Carmen.app/Contents/MacOS/" << cmd.at(0);
            wxString fn = wxString::wxString( ss_fn.str().c_str() );
            if( ! wxFileName(fn).FileExists() ){
                std::stringstream ss_err;
                ss_err << "ERROR running command. Could not execute:\n" << fn.ToStdString() << "\n";
                ss_err << "File does not exist: " << ss.str().c_str();
                wxMessageBox(wxString::FromAscii( ss_err.str().c_str() ), wxString::FromAscii("Error executing command"), wxOK|wxICON_EXCLAMATION, this);
                this->EndModal(wxID_CANCEL);
            }
            if( !( wxFileName( fn ).IsFileExecutable() )){
                wxFileName( fn ).SetPermissions(wxPOSIX_USER_READ | wxPOSIX_USER_WRITE | wxPOSIX_USER_EXECUTE | wxPOSIX_GROUP_READ | wxPOSIX_GROUP_EXECUTE);
            }
        }
        else{
			ss << ".exe";
        }
		ss << "\" "; // trailing space is intentional
	}

    
	for(int i=1; i<(int)cmd.size();i++){
		if( this->investigation->is_mac )
			boost::algorithm::replace_all(cmd.at(i), "\\", "/");
		else
			boost::algorithm::replace_all(cmd.at(i), "/", "\\");
		ss << "\"";
		ss << cmd.at(i);
		ss << "\"";
		if( i<(int)cmd.size()-1){
			ss << " ";
		}
	}
    this->call = ss.str();
    if( investigation->cp->get(std::string("Internal"), std::string("debug")).compare(std::string("true")) == 0 )
    	std::cout << ss.str() << "\n";
    if( this->investigation->is_mac ){
    	// wxMacExecute (utilsexc_base.cpp) is brain-damaged and won't run anything
    	// on the mac unless it's in an application bundle.
    	// This is a work-around using simple unix primitives.
#ifndef WIN32
   		this->pid = popen2(this->call, this->file_descriptor);
#endif
		this->stream = fdopen(this->file_descriptor, "r");
    }
    else{
    	this->process = new wxProcess(this, ID_PROCESS);
    	this->process->Redirect();
    	this->pid = wxExecute( wxString::FromAscii(ss.str().c_str()), wxEXEC_ASYNC, this->process );
    }
    if( this->pid<1){
        std::stringstream ss_err;
        ss_err << "ERROR running command.  Could not execute:\n" << ss.str();
        wxMessageBox(wxString::FromAscii( ss_err.str().c_str() ), wxString::FromAscii("Error executing command"), wxOK|wxICON_EXCLAMATION, this);
        this->EndModal(wxID_CANCEL);
    }
    this->show_progress = show_progress;
    if(show_progress){
    	this->timer = new wxTimer(this, ID_TIMER);
    	timer->Start(10);
    }
    else{
    	CloseDialog();
    }
}


void ProgressDialog::OnClickCancel( wxCommandEvent& event ){
	this->cancelled = true;
	if(this->investigation->is_mac){
#ifndef WIN32
		kill(this->pid, SIGKILL);
#endif
		CloseDialog();
	}
	else{
		if(this->process != NULL)
			this->process->Kill(this->pid, wxSIGKILL);
	}
}

bool ProgressDialog::ReadFromProcess(){

	if(this->investigation->is_mac){
		//wxFileInputStream fileinput(this->file_descriptor);
		if(stream==NULL){
			std::cout.flush();
			wxMessageBox(wxString::FromAscii("Pipe error"));
			return false;
		}
		std::string msg("");
		char line[1000];
		for(int i=0;i<1000;i++)
			line[i]=EOF;

		if( fgets(line, 1000, stream) ){
			msg = std::string(line);
			this->lbx_progress->Append(wxString::FromAscii(msg.c_str()));
			if( msg.substr(0,6).compare("ERROR:")==0)
				this->error = this->error + msg;
		}
		else{
			return false;
		}
	}
	else{
		std::cout.flush();
		while(this->process != NULL && this->process->IsInputAvailable()){
			wxTextInputStream in(*this->process->GetInputStream());
			wxString msg;
			msg << in.ReadLine();
			this->lbx_progress->Append(msg);
			std::string chk(msg.ToAscii());
			if( chk.substr(0,6).compare("ERROR:")==0)
				this->error = this->error + chk;
		}
	}
	int cur_count = this->lbx_progress->GetCount()-1;
	if( cur_count>0 )
		this->lbx_progress->SetFirstItem( cur_count );
	this->Refresh();
	this->Update();
	this->Fit();
	return true;
}

void ProgressDialog::tick(wxTimerEvent& evt){
	if(!ReadFromProcess())
		CloseDialog();
}

void ProgressDialog::OnProcessEnded( wxProcessEvent& WXUNUSED(event) ){
	CloseDialog();
}


void ProgressDialog::CloseDialog(){
	std::string empty("");
	std::string error_pattern("ERROR: ");
	std::string slash_r("\r");
	
	if( !this->cancelled ){
		//ReadFromProcess(false);
		//this->process = NULL;
		if( this->error.size()>0 ){
			std::stringstream errmsg;
			errmsg << "ERROR running command:\n" << this->call << "\n\nReported error:\n" << this->error;
			wxMessageBox(wxString::FromAscii(errmsg.str().c_str()), wxString::FromAscii("Error executing command"), wxOK|wxICON_EXCLAMATION, this);
			this->cancelled = true;
		}
	}
	if(this->timer->IsRunning())
		this->timer->Stop();
	delete this->timer;

    if(this->cancelled){
    	this->EndModal(wxID_CANCEL);
    }
    else{
    	this->EndModal(wxID_OK);
    }
}


void ProgressDialog::OnIdle( wxIdleEvent& evt ){
	if(!ReadFromProcess())
		CloseDialog();
}
