from flask import Flask, render_template, request
from wtforms import Form, StringField, TextAreaField, PasswordField, validators
from Bio.SubsMat import MatrixInfo

app= Flask(__name__)




@app.route("/")
@app.route('/send',methods=['GET','POST'])
def send():
    if request.method == 'POST':
        seq1 = request.form['seq1']
        seq2 = request.form['seq2']
        gap = request.form['gap']
        gap = int(gap)
        seq2Len=len(seq2)
        seq1Len=len(seq1)
        matris = request.form['matris']
        if matris=='blosum62':
            seqmatris = MatrixInfo.blosum62
        elif matris == 'pam120':
            seqmatris = MatrixInfo.pam120
        elif matris == 'pam250':
            seqmatris = MatrixInfo.pam250
        OurMatrix =[[0 for x in range(seq1Len+1)] for y in range(seq2Len+1)]
        for i in range(seq1Len+1):
            if i == 0:
                OurMatrix[0][i]=0
            else:
                OurMatrix[0][i] =-(i)*gap
        for i in range(seq2Len+1):
            if i == 0:
                OurMatrix[i][0] = 0
            else:
                OurMatrix[i][0] = -(i) * gap
        sideMatrix = [[0 for x in range(seq1Len)] for y in range(seq2Len)]
        bestway=[]
        for i in range(len(OurMatrix)):
            for j in range(len((OurMatrix[i]))):
                if i == 0 and j == 0:
                    continue
                if OurMatrix[i][j] == 0:
                    ustdeger=OurMatrix[i-1][j]
                    soldeger = OurMatrix[i][j-1]
                    caprazdeger = OurMatrix[i-1][j-1]
                    ustdeger=ustdeger-gap
                    soldeger=soldeger-gap
                    a=seq1[j-1]
                    b=seq2[i-1]
                    pair = (a,b)
                    if pair not in seqmatris:
                        reverse_pair = tuple(reversed(pair))
                        caprazdeger= seqmatris[reverse_pair]
                    else:
                        caprazdeger=seqmatris[pair]
                    caprazdeger += OurMatrix[i-1][j-1]

                    result = Hesapla(ustdeger,soldeger,caprazdeger)
                    OurMatrix[i][j] = result
                    if result == ustdeger:
                        yon = 'up'
                    if result == soldeger:
                        yon = 'left'
                    if result == caprazdeger:
                        yon = 'cross'
                    sideMatrix[i-1][j-1] = yon



        dizin(seq2Len-1,seq1Len-1,bestway,sideMatrix,seq2)
        dizilim = []
        for i in range(seq1Len):
            dizilim.append(bestway[i])

        return render_template('answer.html', seq1=seq1, seq2=seq2, gap=gap,seqmatris=seqmatris,seq2Len=seq2Len,seq1Len=seq1Len,OurMatrix=OurMatrix,sideMatrix=sideMatrix,dizilim=reversed(dizilim))

    return render_template('form.html')

def Hesapla(a,b,c):
    list1 = [a,b,c];
    return max(list1)

def dizin(i,j,dizi,anaDizi,seq):
    if i>=0 and j>=0:
        if anaDizi[i][j] == 'cross':
            a=seq[i]
            dizi.append(a)
            i-=1
            j-=1
            if i>=0 or j>=0:
                dizin(i, j, dizi, anaDizi, seq)

        if anaDizi[i][j] == 'left':
            dizi.append('-')
            j-=1
            if i>=0 or j>=0:
                dizin(i, j, dizi, anaDizi, seq)


if __name__ == '__main__':
    app.run(debug=True)

