//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
#include <stdio.h>
#include <math.h>
#include "Unit1.h"
#include "Unit2.h"
#include <vcl\Clipbrd.hpp>
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

int start_x = 50, start_y = 50;


void __fastcall TForm1::Button1Click(TObject *Sender)
{
        double IMG_Y[200][500];
        double IMG_Z[200][500];

        FILE *file = fopen("image_profile1.txt", "r");

        long double A = 25.3342L;

        Image1->Width = ceill(A*20.0L);
        Image1->Height = Image1->Width;

        Image1->Canvas->Refresh();

        Image1->Canvas->Pen->Color = clSilver;

        int ii=0;
        while (ii < Image1->Width)
        {
                Image1->Canvas->MoveTo(ii,0);
                Image1->Canvas->LineTo(ii, Image1->Height);

                ii = ii + 20;
        }

        ii=0;
        while (ii < Image1->Height)
        {
                Image1->Canvas->MoveTo(0, ii);
                Image1->Canvas->LineTo(Image1->Width, ii);

                ii = ii + 20;
        }

        for (int i=0; i<200; i++)
        {
                for (int j=0; j<500; j++)
                {
                        fscanf(file, "%le\n", &IMG_Y[i][j]);
                        fscanf(file, "%le\n", &IMG_Z[i][j]);

                        int y = ceill(IMG_Y[i][j]*20.0L);
                        int z = ceill(IMG_Z[i][j]*20.0L);

                        Image1->Canvas->Pixels[ y ][ z ] = clBlack;
                        //Image1->Canvas->Pixels[ y ][ z ] = clRed;
                }
        }

}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button2Click(TObject *Sender)
{
        if ( SavePictureDialog1->Execute() )
        {
                Image1->Picture->SaveToFile( SavePictureDialog1->FileName + ".bmp");
        }
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button3Click(TObject *Sender)
{
        double IMG_Y[200][500];
        double IMG_Z[200][500];

        FILE *file = fopen("image_profile2.txt", "r");

        long double A = 25.3342L;

        Image1->Width = ceill(A*20.0L);
        Image1->Height = Image1->Width;

        Image1->Canvas->Refresh();

        Image1->Canvas->Pen->Color = clSilver;

        int ii=0;
        while (ii < Image1->Width)
        {
                Image1->Canvas->MoveTo(ii,0);
                Image1->Canvas->LineTo(ii, Image1->Height);

                ii = ii + 20;
        }

        ii=0;
        while (ii < Image1->Height)
        {
                Image1->Canvas->MoveTo(0, ii);
                Image1->Canvas->LineTo(Image1->Width, ii);

                ii = ii + 20;
        }


        for (int i=0; i<200; i++)
        {
                for (int j=0; j<500; j++)
                {
                        fscanf(file, "%le\n", &IMG_Y[i][j]);
                        fscanf(file, "%le\n", &IMG_Z[i][j]);

                        int y = ceill(IMG_Y[i][j]*20.0L);
                        int z = ceill(IMG_Z[i][j]*20.0L);

                        Image1->Canvas->Pixels[ y ][ z ] = clBlack;
                        Image1->Canvas->Pixels[ y ][ z ] = clBlue;
                }
        }
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Button4Click(TObject *Sender)
{
        double IMG_Y[200][500];
        double IMG_Z[200][500];

        FILE *file = fopen("image_profile3.txt", "r");

        long double A = 25.3342L;

        Image1->Width = ceill(A*20.0L);
        Image1->Height = Image1->Width;

        Image1->Canvas->Refresh();

        Image1->Canvas->Pen->Color = clSilver;

        int ii=0;
        while (ii < Image1->Width)
        {
                Image1->Canvas->MoveTo(ii,0);
                Image1->Canvas->LineTo(ii, Image1->Height);

                ii = ii + 20;
        }

        ii=0;
        while (ii < Image1->Height)
        {
                Image1->Canvas->MoveTo(0, ii);
                Image1->Canvas->LineTo(Image1->Width, ii);

                ii = ii + 20;
        }


        for (int i=0; i<200; i++)
        {
                for (int j=0; j<500; j++)
                {
                        fscanf(file, "%le\n", &IMG_Y[i][j]);
                        fscanf(file, "%le\n", &IMG_Z[i][j]);

                        int y = ceill(IMG_Y[i][j]*20.0L);
                        int z = ceill(IMG_Z[i][j]*20.0L);

                        Image1->Canvas->Pixels[ y ][ z ] = clBlack;
                        Image1->Canvas->Pixels[ y ][ z ] = clGreen;
                }
        }
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button5Click(TObject *Sender)
{
        double IMG_Y[200][500];
        double IMG_Z[200][500];

        FILE *file = fopen("image_profile4.txt", "r");
        
        for (int i=0; i<200; i++)
        {
                for (int j=0; j<500; j++)
                {
                        fscanf(file, "%le\n", &IMG_Y[i][j]);
                        fscanf(file, "%le\n", &IMG_Z[i][j]);

                        int y = ceill(IMG_Y[i][j]*20.0L);
                        int z = ceill(IMG_Z[i][j]*20.0L);

                        Image1->Canvas->Pixels[ y ][ z ] = clBlack;
                        Image1->Canvas->Pixels[ y ][ z ] = clRed;
                }
        }
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button7Click(TObject *Sender)
{
        Edit1->Text = StringReplace(Edit1->Text, ".", ",", TReplaceFlags() << rfReplaceAll);
        long double A = StrToFloat(Edit1->Text);

        // стираем предыдущее изображение
        Image1->Canvas->Rectangle(-1, -1, Image1->Width + 1, Image1->Width + 1);
        Image1->Width = ceill(A*20.0L) + 2*start_x;
        Image1->Height = Image1->Width;
        Image1->Picture->Bitmap->Height = Image1->Height;
        Image1->Picture->Bitmap->Width = Image1->Width;
        // очищаем область
        Image1->Canvas->Refresh();
        Image1->Canvas->Rectangle(-1, -1, Image1->Width + 1, Image1->Width + 1);
        Image1->Canvas->Pen->Color = clSilver;

        int ii = Image1->Width;
        while (ii > start_x*2)
        {
                Image1->Canvas->MoveTo(start_x + Image1->Width - ii, start_y);
                Image1->Canvas->LineTo(start_x + Image1->Width - ii, Image1->Height - start_y);

                Image1->Canvas->MoveTo(start_x, ii - start_y);
                Image1->Canvas->LineTo(Image1->Width - start_x, ii - start_y);

                ii = ii - 20;
        }

        Image1->Canvas->Pen->Color = clBlack;
        Image1->Canvas->Pen->Width = 2;

        // Рисуем оси координат:
        Image1->Canvas->MoveTo(start_x - 5, Image1->Width - start_y + 2);
        Image1->Canvas->LineTo(Image1->Width - start_x/2, Image1->Width - start_y + 2);
        Image1->Canvas->LineTo(Image1->Width - start_x/2 - 10, Image1->Width - start_y - 3);
        Image1->Canvas->LineTo(Image1->Width - start_x/2, Image1->Width - start_y + 2);
        Image1->Canvas->LineTo(Image1->Width - start_x/2 - 10, Image1->Width - start_y + 7);

        Image1->Canvas->MoveTo(start_x - 1, Image1->Width - start_y + 5);
        Image1->Canvas->LineTo(start_x - 1, start_y/2);
        Image1->Canvas->LineTo(start_x - 6, start_y/2 + 10);
        Image1->Canvas->LineTo(start_x - 1, start_y/2);
        Image1->Canvas->LineTo(start_x + 4, start_y/2 + 10);

        // Делаем подписи к осям координат
        Image1->Canvas->Font->Size = 12;
        Image1->Canvas->TextOutA(start_x - 15, Image1->Width - start_y + 5, "0");
        Image1->Canvas->TextOutA(Image1->Width - start_x/2, Image1->Width - start_y + 5, "Y");
        Image1->Canvas->TextOutA(start_x - 20, start_y/2 - 5, "Z");
        ii = 5;
        while (ii < A+1 ) {
            Image1->Canvas->MoveTo(start_x + 20*ii, Image1->Width - start_y - 2);
            Image1->Canvas->LineTo(start_x + 20*ii, Image1->Width - start_y + 5);

            Image1->Canvas->TextOutA(start_x + 20*ii - 5, Image1->Width - start_y + 5, IntToStr(ii));

            Image1->Canvas->MoveTo(start_x - 5, Image1->Width - start_y - 20*ii);
            Image1->Canvas->LineTo(start_x + 2, Image1->Width - start_y - 20*ii);

            Image1->Canvas->TextOutA(start_x - 25, Image1->Width - start_y - 20*ii - 10, IntToStr(ii));

            ii += 5;
        }

        Image1->Canvas->Pen->Width = 1;

}
//---------------------------------------------------------------------------




void __fastcall TForm1::Button8Click(TObject *Sender)
{
Form2->Show();        
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button6Click(TObject *Sender)
{
        Image1->Picture->Bitmap->Height = Image1->Height;
        Image1->Picture->Bitmap->Width = Image1->Width;
        Clipboard()->Assign(Image1->Picture->Bitmap);     
}
//---------------------------------------------------------------------------

