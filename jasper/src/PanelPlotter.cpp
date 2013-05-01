#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <set>
#include <math.h>
#include "wx/wx.h"
#include "wx/sizer.h"
#include "wx/colour.h"
#include "PanelPlotter.h"
using namespace std;

BEGIN_EVENT_TABLE(Plotter, wxPanel)
EVT_PAINT(Plotter::paintEvent)
END_EVENT_TABLE()

Plotter::Plotter(wxWindow* parent, int id, int w, int h) : wxPanel(parent, id, wxDefaultPosition, wxSize(w,h)), DEFAULT_MARGIN_Y_LOWER(100), DEFAULT_MARGIN_Y_UPPER(50), DEFAULT_MARGIN_X(50) {
	this->parent = parent;
	NA_VALUE = -99999;
	this->panel_width = w;
	this->panel_height = h;
	this->scale_v = 1;
	this->scale_h = 1;
	this->margin_x = this->DEFAULT_MARGIN_X;
	this->margin_y_upper = this->DEFAULT_MARGIN_Y_UPPER;
	this->margin_y_lower = this->DEFAULT_MARGIN_Y_LOWER;
	this->y_orig = this->panel_height - this->margin_y_lower;
	this->pen_width=2;
	this->key_location = KEY_UPPER_RIGHT;
	this->draw_lines_on_plot=true;
	this->colors.push_back( new wxColour(0,0,0) );
	this->colors.push_back( new wxColour(0,120,45) );
	this->colors.push_back( new wxColour(255,0,0) );
	this->colors.push_back( new wxColour(220,0,220) );
	this->colors.push_back( new wxColour(0,150,150) );
	this->colors.push_back( new wxColour(50,0,150) );
	this->colors.push_back( new wxColour(50,220,0) );
	this->colors.push_back( new wxColour(150,150,150) );
	this->colors.push_back( new wxColour(160,0,50) );
	this->colors.push_back( new wxColour(255,150,0) );
}

void Plotter::SetDimensions(int height, int width){
	this->panel_height=height;
	this->y_orig = this->panel_height - this->margin_y_lower;
	this->panel_width=width;
	this->SetSize(width,height);
}

void Plotter::DrawLines(bool dl){
	this->draw_lines_on_plot = dl;
}

void Plotter::paintEvent(wxPaintEvent & evt){
    wxPaintDC dc(this);
    render(dc);
}

void Plotter::Repaint(){
    wxClientDC dc(this);
    render(dc);
}

void Plotter::set_key_location(int keyloc ){
	if( keyloc==KEY_UPPER_RIGHT || keyloc==KEY_UPPER_LEFT || keyloc==KEY_LOWER_LEFT || keyloc==KEY_LOWER_RIGHT)
		this->key_location = keyloc;
	else
		throw std::string("Invalid key location specified");
}

void Plotter::set_scale(){
	// set horizontal scale scale_h, number of pixels between between samples
	// set vertical scale scale_v, number of pixels between two integers
	this->max_y = this->min_y = 0;
	if(this->values.size()==0)
		return;

	double v=0;
	for(int i=0; i<(int)this->values.size(); i++){
		for(int j=0; j<(int)this->values.at(i)->size(); j++){
			v = this->values.at(i)->at(j);
			if(v==NA_VALUE)
				continue;
			if(v>this->max_y)
				this->max_y = v;
			if(v<this->min_y)
				this->min_y = v;
		}
	}
	int nHoriz = (int)this->values.at(0)->size();
	int nVert = (int)ceil(this->max_y) - (int)floor(this->min_y) + 1;
	if(nHoriz<1)
		this->scale_h=0;
	else{
		this->scale_h = ( this->panel_width-(2*this->margin_x) )  / nHoriz;
	}
	if(nVert<1)
		this->scale_v=0;
	else
		this->scale_v = (int) ( ( this->panel_height - (this->margin_y_upper+this->margin_y_lower))  / (nVert) );
}

void Plotter::draw_axis(wxDC& dc){
	dc.SetBrush(*wxBLACK_BRUSH);
	int tic_size = 10;

	wxCoord x1, x2, y1, y2;

	// vertical axis
	x1 = x2 = this->margin_x;
	y1 = this->margin_y_upper;
	y2 = this->y_orig;
	dc.DrawLine(x1, y1, x2, y2);

	// horizontal axis
	// if no Y value is less than zero, draw at the lower left (origin)
	// else, shift up from the origin
	x2 = this->panel_width - this->margin_x;
	y1 = y2;
	int cur_y=0;
	if( this->min_y < 0 ){
		cur_y = (int)this->min_y - 1;
		y1 = y2 = this->y_orig + (cur_y * this->scale_v );
	}
	dc.DrawLine(x1, y1, x2, y2); // horizontal axis

	// if no data, don't draw tic marks
	if( this->values.size() == 0 )
		return;

	// draw vertical tic marks
	x1 = this->margin_x - tic_size;
	x2 = this->margin_x;
	if( this->scale_v > 0 ){
		for(y1=this->y_orig; y1>=this->margin_y_upper; y1=y1-this->scale_v){
			dc.DrawLine(x1, y1, x2, y1);
			std::stringstream ss;
			ss << cur_y++;
			dc.DrawText(wxString::FromAscii(ss.str().c_str()), x1-20, y1-10);
		}
	}
	// draw horizontal tic marks
	y1 = this->y_orig;
	y2 = y1+tic_size;
	x1 = this->margin_x;
	for(int i=0; i<(int)this->sample_labels.size(); i++){
		dc.DrawLine(x1, y1, x1, y2);
		x1 += this->scale_h;
	}
}

void Plotter::reset(){
	this->sample_labels.clear();
	for(int i=0; i<(int)this->values.size(); i++)
		delete values.at(i);
	this->values.clear();
	this->value_labels.clear();
	this->margin_y_lower = this->DEFAULT_MARGIN_Y_LOWER; // can be altered if sample names are long
}


void Plotter::add_valueset(std::string name, std::vector<double>* values_in){
	if( this->values.size()>0 ){
		int cur_n = (int)values_in->size();
		if( (int)values.at(0)->size() != cur_n ){
			throw std::string("New valueset length does not match first valueset length");
		}
	}
	this->value_labels.push_back(name);
	std::vector<double>* V = new std::vector<double>();
	for(int i=0; i<(int)values_in->size(); i++)
		V->push_back(values_in->at(i));
	this->values.push_back(V);
}

void Plotter::set_plot_values(std::vector<std::string> sample_labels, std::vector<std::string> gene_labels, std::vector<std::vector<double>*>& values){
	this->reset();
	if(gene_labels.size() != values.size())
		throw std::string("gene_labels and values do not have compatible size");
	for(int i=0; i<(int)gene_labels.size(); i++){
		if(sample_labels.size() != values.at(i)->size() )
			throw std::string("sample labels and values do not have compatible size");
		add_valueset( gene_labels.at(i), values.at(i));
	}
	int max_label_length=0;
	for(int i=0; i<(int)sample_labels.size(); i++){
		this->sample_labels.push_back(sample_labels.at(i));
		if( (int)sample_labels.at(i).size() > max_label_length )
			max_label_length = (int)sample_labels.at(i).size();
	}
	if( max_label_length > 14 ){
		this->margin_y_lower += 50;
		this->y_orig = this->panel_height - this->margin_y_lower;
	}
}


void Plotter::draw_plot(wxDC& dc){
	int x_cur;
	std::vector<double>* V;
	for(int i=0;i< (int)values.size(); i++){
		V = values.at(i);
		dc.SetBrush(*this->colors.at( i % (int)this->colors.size() ));
		wxPen pen(*this->colors.at( i % (int)this->colors.size() ));
		pen.SetWidth(this->pen_width);
		if(i>=(int)this->colors.size())
			pen.SetStyle(wxLONG_DASH);
		dc.SetPen(pen);
		x_cur = this->margin_x;
		if( (int)V->size()==1){
			// special case of only one sample; next block assumes at least two samples.
			wxCoord y = this->y_orig - (V->at(0) * this->scale_v);
			if(this->min_y < 0 ){
				int offset = (int)floor(this->min_y) - 1;
				y += (offset * this->scale_v );
			}
			if( V->at(0)!=NA_VALUE )
				dc.DrawCircle(x_cur, y, 2);
		}
		else{
			for(int j=0; j<(int)V->size()-1; j++){
				wxCoord y1 = this->y_orig - (V->at(j) * this->scale_v);
				wxCoord y2 = this->y_orig - (V->at(j+1) * this->scale_v);
				if(this->min_y < 0 ){
					int offset = (int)this->min_y - 1;
					y1 += (offset * this->scale_v );
					y2 += (offset * this->scale_v );
				}
				wxCoord x1 = x_cur;
				wxCoord x2 = x1 + this->scale_h;
				if( this->draw_lines_on_plot ){
					if( V->at(j)!=NA_VALUE && V->at(j+1)!= NA_VALUE)
						dc.DrawLine(x1, y1, x2, y2);
				}
				if( V->at(j)!=NA_VALUE )
					dc.DrawCircle(x1, y1, 2);
				if( V->at(j+1)!=NA_VALUE )
					dc.DrawCircle(x2, y2, 2);
				x_cur += this->scale_h;
			}
		}
	}
}

void Plotter::draw_sample_labels(wxDC& dc){
	int x_cur = this->margin_x;
	dc.SetFont(wxFont(8,wxFONTFAMILY_SWISS,wxFONTSTYLE_NORMAL,wxFONTWEIGHT_BOLD));
	for(int i=0; i<(int)this->sample_labels.size(); i++){
		int nchar = this->sample_labels.at(i).size();
		int y_cur = this->y_orig;
		if( nchar < 8 )
			y_cur = y_cur+50;
		else if( nchar < 20 )
			y_cur += 80;
		else
			y_cur += 95;

		dc.DrawRotatedText( wxString::FromAscii(this->sample_labels.at(i).c_str()), x_cur-5, y_cur, 90);
		x_cur += this->scale_h;
	}
}

void Plotter::draw_key(wxDC& dc){

	// if no data, don't draw tic marks
	if( this->values.size() == 0 )
		return;

	int x0 = this->panel_width - this->margin_x - 100;
	int y0 = this->margin_y_upper + 10;

	if( this->key_location==KEY_UPPER_RIGHT){
		x0 = this->panel_width - this->margin_x - 100;
		y0 = this->margin_y_upper + 10;
	}
	else if(this->key_location==KEY_UPPER_LEFT){
		x0 = this->margin_x + 50;
		y0 = this->margin_y_upper + 10;
	}
	else if(this->key_location==KEY_LOWER_RIGHT){
		x0 = this->panel_width - this->margin_x - 100;
		y0 = this->panel_height - this->margin_y_lower - 100;
	}
	else if(this->key_location==KEY_LOWER_LEFT){
		x0 = this->margin_x + 50;
		y0 = this->panel_height - this->margin_y_lower - 100;
	}

	int y_cur = y0;
	dc.SetFont(wxFont(8,wxFONTFAMILY_SWISS,wxFONTSTYLE_NORMAL,wxFONTWEIGHT_BOLD));
	for(int i=0; i<(int)this->value_labels.size(); i++){
		dc.DrawText( wxString::FromAscii(this->value_labels.at(i).c_str()), x0, y_cur);
		wxPen pen(*this->colors.at( i % (int)this->colors.size() ));
		pen.SetWidth(6);
		dc.SetPen(pen);
		if(i>=(int)this->colors.size()){
			pen.SetWidth(3);
			dc.DrawLine(x0-15, y_cur+5, x0-14, y_cur+5);
			dc.DrawLine(x0-7,y_cur+5, x0-5, y_cur+5);
		}
		else
			dc.DrawLine(x0-15, y_cur+5, x0-5, y_cur+5);
		y_cur += 10;
	}
	wxPen pen(*wxBLACK_PEN);
	dc.SetBrush(*wxTRANSPARENT_BRUSH);
	dc.SetPen(pen);
	dc.DrawRectangle(x0-25, y0-5, 150, y_cur - y0 + 10 );
}

void Plotter::render(wxDC&  dc){
	dc.SetBackground( *wxWHITE_BRUSH );
	set_scale();
	draw_axis(dc);
	draw_plot(dc);
	draw_sample_labels(dc);
	draw_key(dc);
}

Plotter::~Plotter(){
	this->reset();
}
