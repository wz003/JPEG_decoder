//
//  main.c
//  JPEG_decoder
//
//  Created by 至 on 2017/4/30.
//  Copyright © 2017年 至. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>



///
const char hex_chars[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };

int corder[64] = {0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,4,5,6,7,7,6,5,4,3,2,1,2,3,4,5,6,7,7,6,5,4,3,4,5,6,7,7,6,5,6,7,7};
int rorder[64] = {0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,7,6,5,4,3,2,3,4,5,6,7,7,6,5,4,5,6,7,7,6,7};

struct node {
    node* left;
    node* right;
    int data;
};
int width=0;
int height=0;
int y_sample[2] = {0};
node *tree_DC_00 = new node;
node *tree_AC_00 = new node;
node *tree_DC_01 = new node;
node *tree_AC_01 = new node;
int q_table[2][8][8] = {0};


void IDCT2(double *F, double *f)
{
    const double pi = acos(-1);
    const double c30 = (2*cos(pi/8));
    const double c45 = sqrt(2);
    const double c60 = (2*sin(pi/8));
    const double Q = c30 - c60;
    const double R = c30 + c60;
    double tmp1;
    double m[8];
    
    m[0] = F[0];
    m[1] = F[4];
    m[2] = F[2] - F[6];
    m[3] = F[2] + F[6];
    m[4] = F[5] - F[3];
    m[5] = F[1] + F[7] - F[3] - F[5];
    m[6] = F[1] - F[7];
    m[7] = F[1] + F[7] + F[3] + F[5];
    
    tmp1 = c60 * (m[4] + m[6]);
    F[0] = m[0];
    F[1] = m[1];
    F[2] = m[2] * c45;
    F[3] = m[3];
    F[4] = -Q * m[4] - tmp1;
    F[5] = m[5] * c45;
    F[6] = R * m[6] - tmp1;
    F[7] = m[7];
    
    m[0] = F[6] - F[7] - F[5];
    m[1] = F[0] - F[1];
    m[2] = F[2] - F[3];
    m[3] = F[0] + F[1];
    m[4] = F[6] - F[7];
    m[5] = F[4];
    m[6] = F[3];
    m[7] = F[7];
    
    F[0] = m[7];
    F[1] = m[0];
    F[2] = m[4];
    F[3] = m[1] + m[2];
    F[4] = m[3] + m[6];
    F[5] = m[1] - m[2];
    F[6] = m[3] - m[6];
    F[7] = m[5] - m[0];
    
    f[0*8] = F[4] + F[0];
    f[1*8] = F[3] + F[2];
    f[2*8] = F[5] - F[1];
    f[3*8] = F[6] - F[7];
    f[4*8] = F[6] + F[7];
    f[5*8] = F[5] + F[1];
    f[6*8] = F[3] - F[2];
    f[7*8] = F[4] - F[0];
}

void IDCT(int a[8][8])
{
    const double pi = acos(-1);
    const double X[8] =
    {
        sqrt(2)*cos(0 *pi/16)/4,
        cos(1*pi/16)/2,
        cos(2*pi/16)/2,
        cos(3*pi/16)/2,
        cos(4*pi/16)/2,
        cos(5*pi/16)/2,
        cos(6*pi/16)/2,
        cos(7*pi/16)/2
    };
    
    int i, j;
    double x[64], y[64];
    
    for(i = 0; i < 8; i++)
        for(j = 0; j < 8; j++)
            x[i*8+j] = a[i][j] * X[i] * X[j];
    
    for(i = 0; i < 8; i++)
        IDCT2(x+8*i, y+i);
    
    for(i = 0; i < 8; i++)
        IDCT2(y+8*i, x+i);
    
    for(i = 0; i < 8; i++)
        for(j = 0; j < 8; j++)
            a[i][j] = (int)(x[i*8+j] + 0.5);
}


void Idct(int in[8][8], int out[8][8])
{
    int i, j, u, v;
    double s;
    
    for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++)
        {
            s = 0;
            
            for (u = 0; u < 8; u++)
                for (v = 0; v < 8; v++)
                    s += in[u][v] * cos((2 * i + 1) * u * M_PI / 16) *
                    cos((2 * j + 1) * v * M_PI / 16) *
                    ((u == 0) ? 1 / sqrt(2) : 1.) *
                    ((v == 0) ? 1 / sqrt(2) : 1.);
            
            out[i][j] =round (s / 4)+128;
        }
}

void Quna_Table(unsigned char *buffer , int i)
{
    //buffer[i,i+1] : FFDB
    i+=2;
    //int len = buffer[i]*256+buffer[i+1];
    //printf("Quan_Table length is : %d\n",len);
    i+=2;
    
    //int precision = buffer[i]/16; //0>>8bits , 1>>16bits
    //int id = buffer[i]%16; // 量化表id(0~3)
    //printf("precision:%dbits,id:%d\n",(precision+1)*8,id);
    i++;
    
    int count = 0;
    while (count <2) {
        int temp[64] = {0};
        for (int j = 0;j<64;j++) {
            temp[j] = buffer[i];
            i++;
        }
        
        for (int j = 0;j<64;j++) {
            q_table[count][corder[j]][rorder[j]] = temp[j];
        }
        if (buffer[i] == 255) {
            i+=4;
        }
        i++;
        count++;
    }
//    for (int l = 0;l<2;l++) {
//        printf("%dth quan table : \n",l);
//        for (int j = 0;j<8;j++) {
//            for (int k = 0;k<8;k++) {
//                printf("%d ",q_table[l][j][k]);
//            }
//            printf("\n");
//        }
//        printf("\n");
//        
//    }
    
}



void SOF(unsigned char *buffer , int i)
{
    height = buffer[i+5]*256+buffer[i+6];
    width = buffer[i+7]*256+buffer[i+8];
    y_sample[0] = buffer[i+11]/16;
    y_sample[1] = buffer[i+11]%16;
//    printf("resolution:height:%d,width:%d\n",height,width);
//    printf("y_sample:%d,%d\n",y_sample[0],y_sample[1]);
}

void DHT(unsigned char *buffer , int i)
{
    
    
    int t_num[2][2][16]={0};
    
    int count = 0;
    int len;
    //FFC4 = 65476
    int sp = 0;
    
    
    while (buffer[i]*256+buffer[i+1] != 65498) {//還沒到SOC(FFDA)區間
        
        if ( buffer[i]*256+buffer[i+1] == 65476 ||sp == 1) {
            //buffer[i],buffer[i+1] : FFC4
            //printf("Huffman_Tree%d:\n",count);
            
            if (sp == 0) {//沒有標記的圖中不會有長度,直接從codeword length開始
                i = i+2;
            
                //printf("strat i is : %d\n",i);
                //buffer[i],buffer[i+1] : length 001F
                len = buffer[i]*256+buffer[i+1];
                //printf("length is : %d\n",len);
                // DC or AC , table ID
                
                
                i = i+2;
                //buffer[i] : DC/AC ,TABLE ID
            }
            
            int type = buffer[i]/16;
            int ID = buffer[i]%16;
            //printf("type:%d,id:%d\n",type,ID);
            //buffer[i] : huffman codeword length counter
            for (int j = 0;j<16;j++) {
                i++;
                t_num[type][ID][j] = buffer[i];
                //printf("%d:%d\n",j+1,t_num[type][ID][j]);
            }
            
            i++;
            //buffer[i] : huffman codeword length counter
            
            
    //            while (i < start+len+2) {
    //                printf("%02X\n",buffer[i]);
    //                i++;
    //            }
            
            //開始造樹
            node *root = new node;
            root->left = NULL;
            root->right = NULL;
            root->data = -1;
            node *current = root;
            
            
            char codeword[] = "0000000000000000";
            
            
            int len_last=0;
            int dif=0;
            int flag = 0; //flag==0 >>長度沒有變動
            bool first = true;
            //len_last : code length (j+1)
            for (int j = 0;j<16;j++) {
                flag = 1;
                //codeword長度變動
                while (t_num[type][ID][j] > 0 ) {
                    
                    //造codeword
                    if (flag == 1) {
                        if (not first) {
                            //只有第一次不用先+1
                            for (int l = 0;l<16;l++) {
                                if (codeword[15-l] == '1' ) {
                                    codeword[15-l] = '0' ;
                                }else {
                                    codeword[15-l] = '1' ;
                                    break;
                                }
                            }
                        }
                        //codeword shift left
                        dif = j+1 - len_last;
                        for (int l =dif;l<16;l++) {
                            codeword[l-dif] = codeword[l];
                        }
                        for (int l = 15;l>15-dif;l--) {
                            codeword[l] = '0' ;
                        }
                        flag = 0;
                        first = false;
                    }else {//flag == 0
                        //codeword length不變 , codeword++
                        for (int l = 0;l<16;l++) {
                            if (codeword[15-l] == '1' ) {
                                codeword[15-l] = '0' ;
                            }else {
                                codeword[15-l] = '1' ;
                                break;
                            }
                        }
                    }
                    
                    len_last = j+1;//記錄現在codeword length
                    
    //                    for (int k=16-len_last;k<16;k++) {
    //                        printf("%c",codeword[k]);
    //                    }
    //                    printf("\n");
                    //codeword 製造完成
                    
                    current = root;
                    
                    for (int k=16-len_last;k<16;k++) {
                        //codeword 為codeword[ (15-len_last+1 ) ~ 15 ] (長為codewordlenght)
                        if (strncmp(&codeword[k], "0", 1)==0 ) {
                            //str = '0' >> 往左走
                            if ( current->left ==NULL) {
                                node *newnode = new node;
                                current->left = newnode;
                                newnode->left = NULL;
                                newnode->right = NULL;
                                newnode->data = -1;
                                current = newnode;
                            }else {
                                current = current->left;
                            }
                        }else {
                            //str = '1' >> 往右走
                            if (current->right == NULL) {
                                node *newnode = new node;
                                current->right = newnode;
                                newnode->left = NULL;
                                newnode->right = NULL;
                                newnode->data = -1;
                                current = newnode;
                            }else {
                                current = current->right;
                            }
                        }
                        //printf("%c",codeword[k]);
                    }
                    //printf(" value: 0x%02x \n",buffer[i]);
                    //這兩個print可以檢查tree structure
                    
                    //已經沿著codeword走到底了
                    
                    
                    current->data = buffer[i];
                    i++;
                    t_num[type][ID][j]--;//造完一個node
                }//這個長度的codeword都處理完了
            }//for迴圈中止,完成一棵tree
            if (buffer[i]*256+buffer[i+1] == 65476) {
                sp = 0;//下一棵樹有FFDA標記
            }else {
                sp = 1;//下一棵樹沒有
            }
            //printf("tree end , next i is :%d\n",i);
            if (type == 0 & ID == 0) {
                tree_DC_00 = root;
            } else if (type ==1 & ID == 0) {
                tree_AC_00 = root;
            } else if (type == 0 & ID ==1) {
                tree_DC_01 = root;
            } else if (type == 1 & ID ==1) {
                tree_AC_01 = root;
            }
                
            
            count = count+1;//count : tree number
            
            
        }else {//沒有讀到tree mark , 繼續往下找
            i = i+1;
        }
    }//找完所有的tree了
}

int MCU_block2[2400][8][16][3] = {0};
void SOS(unsigned char *buffer , int i)
{
    //buffer[i],buffer[i+1] : FFDA
    i+=14;
    //每一張圖這14 bytes皆相同 懶得parse
    //buffer[i]開始為data stream
    int mcu_block[8][8] = {0};
    int used_num = 8;//buffer中已用幾個元素
    int bits[8] = {0};
    int DCorAD = 0;
    int DC_last = 0;
    int y_DC_last = 0;
    int cb_DC_last = 0;
    int cr_DC_last = 0;
    int mcuy[64] = {0};
    int mcu_counter = 0;
    int done[8][8] = {0};
    node *tree = new node;
    node *current = new node;
 
    int h_val = 0;
    int l_val = 0;
    
    //int w = width;
    //int h = height;
    
    int unit_num = 6;
    if (y_sample[1]==1) {
        unit_num = 4;
    }
    int unit_count = 0;
    int unit[unit_num][8][8];
    int MCU_block[1200][16][16][3] = {0};
    //int MCU_block2[2400][8][16][3] = {0};
    int img[480][640][3] = {0}; // 不能用w,h宣告矩陣＠＿＠？
    int total=0;
    while (buffer[i]*256+buffer[i+1] != 65497) { //將所有byte讀出
        
        while (mcu_counter <64) {
            int stream[100] = {0};
            if (DCorAD == 0) { //還沒取出DC
                //printf("unit_count :%d %d\n" , unit_count/unit_num,unit_count%unit_num);
                if ( (unit_count%unit_num)<unit_num-2) {
                    DC_last = y_DC_last;
                    tree = tree_DC_00;
                    current = tree_DC_00;
                }else {
                    if ((unit_count%unit_num)==unit_num-2){
                        DC_last=cr_DC_last;
                    }else if ((unit_count%unit_num)==unit_num-1){
                        DC_last=cb_DC_last;
                        
                    }
                    //printf("DC last is : %d\n",DC_last);
                    tree = tree_DC_01;
                    current = tree_DC_01;
                }
            
                while (current->data == -1) {
                    //開始查tree
                    if (used_num==8) {
                        
                        for (int j = 0;j<8;j++) {
                            bits[j] = buffer[i]>>(7-j) &0x01;
                        }
                        used_num = 0;
                        if (buffer[i]*256+buffer[i+1] != 65280) {
                            i++;
                        }else {
                            i+=2;
                        }
                    }
                    //若buffer[i]使用完,則讀buffer[i+1]
                    
                    if (bits[used_num] == 0) {
                        current = current->left;
                        used_num++;
                    }else {
                        current = current->right;
                        used_num++;
                    }
                    
                    //讀出1bits , 往下走
                    
                }//讀到node data 不為-1
                int stream_len = current->data;
                //printf("stream length : %d\n",stream_len);
                
                for (int j = 0;j<stream_len;j++) {
                    
                    if (used_num==8) {
                        for (int j = 0;j<8;j++) {
                            bits[j] = buffer[i]>>(7-j) &0x01;
                        }
                        used_num = 0;
                        if (buffer[i]*256+buffer[i+1] != 65280) {
                            i++;
                        }else {
                            i+=2;
                        }
                    }
                    //若buffer[i]使用完,則讀buffer[i+1]
                    
                    stream[j] = bits[used_num];
                    used_num++;
                    
                }//往後讀stream_len個bits , 得到stream
                
                int val = 0;
                
                if (stream[0]==0) {
                    for (int j= 0;j<stream_len;j++) {
                        if (stream[j]==0) {
                            stream[j]=1;
                        }else {
                            stream[j] = 0;
                        }
                    }
                    for (int j= 0;j<stream_len;j++) {
                        val += stream[j]<<(stream_len-1-j);
                    }
                    val = -val;
                }else {
                    for (int j= 0;j<stream_len;j++) {
                        val += stream[j]<<(stream_len-1-j);
                    }
                }
                DCorAD = 1;
                
                mcuy[mcu_counter] = val+DC_last;
                
                if ((unit_count%unit_num)<unit_num-2) {
                    y_DC_last= mcuy[mcu_counter];
                    //printf("y_DC_last:%d\n",y_DC_last);
                }else if ((unit_count%unit_num)==unit_num-2) {
                    cr_DC_last=mcuy[mcu_counter];
                }else if ((unit_count%unit_num)==unit_num-1){
                    cb_DC_last=mcuy[mcu_counter];
                }
                
                mcu_counter ++;
                
//                printf("val : %d\n",val);
//                printf("LAST DC : %d\n",DC_last);
                
                //讀出DC值
            }else { // 接著讀出63個AC值
                if ( unit_count%unit_num<unit_num-2) {
                    tree = tree_AC_00;
                    current = tree_AC_00;
                }else {
                    tree = tree_AC_01;
                    current = tree_AC_01;
                }
                while (current->data == -1) {
                    //開始查tree
                    if (used_num==8) {
                        for (int j = 0;j<8;j++) {
                            bits[j] = buffer[i]>>(7-j) &0x01;
                        }
                        used_num = 0;
                        if (buffer[i]*256+buffer[i+1] != 65280) {
                            i++;
                        }else {
                            i+=2;
                        }
                    }
                    //若buffer[i]使用完,則讀buffer[i+1]
                    
                    if (bits[used_num] == 0) {
                        current = current->left;
                        used_num++;
                    }else {
                        current = current->right;
                        used_num++;
                    }
                    
                    //讀出1bits , 往下走
                    
                }//讀到node data 不為-1
                int stream_len = current->data;
                if(stream_len == 0) {//MCU 後面皆為0了
                    for(int j = mcu_counter;j<64;j++) {
                        mcuy[j] = 0;
                    }
                    mcu_counter = 64;
                }else {
                    h_val = stream_len/16;
                    l_val = stream_len%16;
                
                    stream_len = l_val;
                    
                    for (int j = 0;j<stream_len;j++) {
                        
                        if (used_num==8) {
                            for (int j = 0;j<8;j++) {
                                bits[j] = buffer[i]>>(7-j) &0x01;
                            }
                            used_num = 0;
                            if (buffer[i]*256+buffer[i+1] != 65280) {
                                i++;
                            }else {
                                i+=2;
                            }
                        }
                        //若buffer[i]使用完,則讀buffer[i+1]
                        
                        stream[j] = bits[used_num];
                        used_num++;
                        
                    }//往後讀stream_len個bits , 得到stream

                    int val = 0;
                
                    if (stream[0]==0) {
                        for (int j= 0;j<stream_len;j++) {
                            if (stream[j]==0) {
                                stream[j]=1;
                            }else {
                                stream[j] = 0;
                            }
                        }
                        for (int j= 0;j<stream_len;j++) {
                            val += stream[j]<<(stream_len-1-j);
                        }
                        val = -val;
                    }else {
                        for (int j= 0;j<stream_len;j++) {
                            val += stream[j]<<(stream_len-1-j);
                        }
                    }
                    //printf("%dth AC pair is %d,%d\n",mcu_counter,h_val,val);
                    for (int j = 0;j<h_val;j++) {
                        mcuy[mcu_counter] = 0;
                        mcu_counter++;
                    }
                    mcuy[mcu_counter] = val;
                    mcu_counter++;
                }//讀出一對AC值
            
            }
            
        }//讀完一個8*8block
        
        
        
        if (mcu_counter ==64) {
            
            //anti-zig-zag
            for (int j = 0;j<64;j++) {
                mcu_block[corder[j]][rorder[j]] = mcuy[j];
            }
            //de-quanti
            int block[8][8]={0};
            for (int j=0;j<8;j++) {
                for (int k=0;k<8;k++) {
                    if (unit_count%unit_num <unit_num-2) {
                         block[j][k] = mcu_block[j][k] * q_table[0][j][k];
                    }else {
                        block[j][k] = mcu_block[j][k] * q_table[1][j][k];
                    }
                }
            }
            
            //IDCT
            
            //Idct(block, done);
            
            IDCT(block);
            for (int j=0;j<8;j++){
                for (int k=0;k<8;k++){
                    done[j][k]=block[j][k]+128;
                }
            }
            
            //printf("unit_count : %d\n",unit_count%unit_num);
            for (int j=0;j<8;j++) {
                for (int k=0;k<8;k++) {
                    unit[unit_count%unit_num][j][k] = done[j][k];
                    //printf("%4d",unit[unit_count%unit_num][j][k]);
                }
                //printf("\n");
            }
//
            DCorAD = 0;
            mcu_counter = 0;
            
            unit_count++;
            //printf("unit_count:%d\n",unit_count);
            i--;
        }
        
        
        
        
        if  ( unit_count!=0 & unit_count%unit_num==0) {//集滿一組
            //printf("unit_count:%d\n",unit_count);
            total = ((unit_count-1)/unit_num);
            //printf("total : %d\n",total);
            int block16[16][16][3] = {0};
                    //y分量
                if (unit_num ==6) {
                    for(int m = 0;m<8;m++){
                        for (int n =0;n<8;n++) {
                            block16[m][n][0] = unit[0][m][n];
                            block16[m][n+8][0] = unit[1][m][n];
                            block16[m+8][n][0] = unit[2][m][n];
                            block16[m+8][n+8][0] = unit[3][m][n];
                            
                            block16[2*m][2*n][1] = unit[4][m][n];
                            block16[2*m][2*n+1][1] = unit[4][m][n];
                            block16[2*m+1][2*n][1] = unit[4][m][n];
                            block16[2*m+1][2*n+1][1] = unit[4][m][n];
                            
                            block16[2*m][2*n][2] = unit[5][m][n];
                            block16[2*m][2*n+1][2] = unit[5][m][n];
                            block16[2*m+1][2*n][2] = unit[5][m][n];
                            block16[2*m+1][2*n+1][2] = unit[5][m][n];
                            
                        }
                        
                    }
                    
                    for (int k =0;k<16;k++) {
                        for (int j=0;j<16;j++) {
                            for (int l =0;l<3;l++)
                                MCU_block[total][k][j][l] = block16[k][j][l];
                        }
                    }
                }else {//4.jpg
                    
                    for(int m = 0;m<8;m++){
                        for (int n =0;n<8;n++) {
                            block16[m][n][0] = unit[0][m][n];
                            block16[m][n+8][0] = unit[1][m][n];
                         
                            block16[m][2*n][1] = unit[2][m][n];
                            block16[m][2*n+1][1] = unit[2][m][n];
                            
                            block16[m][2*n][2] = unit[3][m][n];
                            block16[m][2*n+1][2] = unit[3][m][n];
                            
                        }
                    }
                    for (int k =0;k<8;k++) {
                        for (int j=0;j<16;j++) {
                            for (int l =0;l<3;l++){
                                
                                MCU_block2[total][k][j][l] = block16[k][j][l];
                            }
                        }
                    }
                }

        }
        
        i++;
    }//讀完所有SOS資料
    

    int w_b =(width/16);
    if (width%16>0) {
        w_b++;
    }
    if (unit_num ==6) {
        for (int j=0;j<total+1;j++) {
            int hor = j%w_b;
            int ver = j/w_b;
            //printf("now:%d,hor:%d,ver:%d\n",j,hor,ver);
            for (int k=0;k<16;k++) {
                for (int l=0;l<16;l++) {
                    img[ver*16+k][hor*16+l][0] = MCU_block[j][k][l][0];
                    img[ver*16+k][hor*16+l][1] = MCU_block[j][k][l][1];
                    img[ver*16+k][hor*16+l][2] = MCU_block[j][k][l][2];
                }
            }
        }
    }else {
        for (int j=0;j<total+1;j++) {
            int hor = j%w_b;
            int ver = j/w_b;
            //printf("now:%d,hor:%d,ver:%d\n",j,hor,ver);
            for (int k=0;k<8;k++) {
                for (int l=0;l<16;l++) {
                    img[ver*8+k][hor*16+l][0] = MCU_block2[j][k][l][0];
                    img[ver*8+k][hor*16+l][1] = MCU_block2[j][k][l][1];
                    img[ver*8+k][hor*16+l][2] = MCU_block2[j][k][l][2];
                }
            }
        }
    }
    
    int temp_RGB[3];
    for (int j=0;j<height;j++) {
        for (int k=0;k<width;k++) {
            
            temp_RGB[0] = img[j][k][0]+ 1.402*(img[j][k][2]-128); //R
            temp_RGB[1] = img[j][k][0]- 0.34414*(img[j][k][1]-128) - 0.71414*(img[j][k][2]-128) ;//G
            temp_RGB[2] = img[j][k][0]+ 1.772*(img[j][k][1]-128); //B
            for (int l = 0;l<3;l++) {
                if (temp_RGB[l]>255) {
                    temp_RGB[l]=255;
                }else if (temp_RGB[l]<0){
                    temp_RGB[l]=0;
                }
                img[j][k][l] = temp_RGB[l];
            }
        }
    }
    
    if (width==286) {
        width = 284;
    }
    
    typedef struct                       /**** BMP file header structure ****/
    {
        unsigned int   bfSize;           /* Size of file */
        unsigned short bfReserved1;      /* Reserved */
        unsigned short bfReserved2;      /* ... */
        unsigned int   bfOffBits;        /* Offset to bitmap data */
    } BITMAPFILEHEADER;
    
    typedef struct                       /**** BMP file info structure ****/
    {
        unsigned int   biSize;           /* Size of info header */
        int            biWidth;          /* Width of image */
        int            biHeight;         /* Height of image */
        unsigned short biPlanes;         /* Number of color planes */
        unsigned short biBitCount;       /* Number of bits per pixel */
        unsigned int   biCompression;    /* Type of compression to use */
        unsigned int   biSizeImage;      /* Size of image data */
        int            biXPelsPerMeter;  /* X pixels per meter */
        int            biYPelsPerMeter;  /* Y pixels per meter */
        unsigned int   biClrUsed;        /* Number of colors used */
        unsigned int   biClrImportant;   /* Number of important colors */
    } BITMAPINFOHEADER;
    
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    
    /* Magic number for file. It does not fit in the header structure due to alignment requirements, so put it outside */
    unsigned short bfType=0x4d42;
    bfh.bfReserved1 = 0;
    bfh.bfReserved2 = 0;
    bfh.bfSize = 2+sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)+width*height;
    bfh.bfOffBits = 0x36;
    
    bih.biSize = sizeof(BITMAPINFOHEADER);
    bih.biWidth = width;
    bih.biHeight = height;
    bih.biPlanes = 1;
    bih.biBitCount = 24;
    bih.biCompression = 0;
    bih.biSizeImage = 0;
    bih.biXPelsPerMeter = 5000;
    bih.biYPelsPerMeter = 5000;
    bih.biClrUsed = 0;
    bih.biClrImportant = 0;
    
    FILE *file = fopen("output.bmp", "wb");
    if (!file)
    {
        printf("Could not write file\n");
        return;
    }
    
    /*Write headers*/
    fwrite(&bfType,1,sizeof(bfType),file);
    fwrite(&bfh, 1, sizeof(bfh), file);
    fwrite(&bih, 1, sizeof(bih), file);
    
    /*Write bitmap*/
    for (int y = bih.biHeight-1; y>=0; y--) /*Scanline loop backwards*/
    {
        for (int x = 0; x < bih.biWidth; x++) /*Column loop forwards*/
        {
            /*compute some pixel values*/
            unsigned char r = img[y][x][0];
            unsigned char g = img[y][x][1];
            unsigned char b = img[y][x][2];
            fwrite(&b, 1, 1, file);
            fwrite(&g, 1, 1, file);
            fwrite(&r, 1, 1, file);
        }
    }
    fclose(file);

    
}



int main(int argc, const char * argv[]) {
    // insert code here...
    FILE *jpg_img;
    long sz;
    
    jpg_img = fopen("2.jpg","rb");
    if (jpg_img == NULL) {
        printf("CANNOT OPEN IMG\n");
        return 1;
    } else {
        printf("File read successfully\n");
        fseek(jpg_img, 0L, SEEK_END);
        sz = ftell(jpg_img);
        
    }
    
    unsigned char buffer[sz];
    
    fseek(jpg_img, 0, SEEK_SET);
    fread(buffer,sizeof(buffer), 1, jpg_img);
    
    int i = 0;
    int x[4];
    char temp[4];
    int part = 0;
    
    
    while( i <sizeof(buffer)) {
        x[0] = buffer[i]/16;
        x[1] = buffer[i]%16;
        x[2] = buffer[i+1]/16;
        x[3] = buffer[i+1]%16;
        temp[0] = hex_chars[x[0]];
        temp[1] = hex_chars[x[1]];
        temp[2] = hex_chars[x[2]];
        temp[3] = hex_chars[x[3]];
        char str[] ={temp[0],temp[1],temp[2],temp[3]};
        
        
        
        if (part == 0 & strncmp(str, "FFD8", 4) ==0 ) {
            printf("SOI，Start of Image at %d\n",i);
            part=1;
        }
        if (part <= 1 & strncmp(str, "FFE0", 4) ==0 ) {
            printf("APP0，Application at %d\n",i);
            
            part=2;
        }
        
        if (part <= 2 & strncmp(str, "FFDB", 4) ==0 ) {
            printf("Define Quantization Table at %d\n",i);
            Quna_Table(buffer,i);
            //check ok;
            part=3;
        }
        
        if (part <= 3 & strncmp(str, "FFC0", 4) ==0 ) {
            printf("Start of Frame at %d\n",i);
            //ok
            
            SOF(buffer , i );
            // resolution
            
            part=4;
        }
        if (part <= 4 & strncmp(str, "FFC4", 4) ==0 ) {
            printf("Define Huffman Table at %d\n",i);
            
            part = 0;
            DHT(buffer , i);
            //check ok

        }
        if (part <= 5 & strncmp(str, "FFDD", 4) ==0 ) {
            printf("Define Restart Interval at %d\n",i);
            part=6;
            
        }
        //ok
        if (part <= 6 & strncmp(str, "FFDA", 4) ==0 ) {
            printf("Start of Scan at %d\n",i);
            part=7;
//
            SOS(buffer , i);
            
        }
        if (part <= 7 & strncmp(str, "FFD9", 4) ==0 ) {
            printf("End of Image at %d\n",i);
            

        }
        
        
        i++;
    }
//    }
    
    fclose(jpg_img);
    
    //hello();
    
    return 0;
}
