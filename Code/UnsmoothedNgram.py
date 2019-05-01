import re
import collections
import unicodedata

from collections import Counter
import math


singleCount={}

def convertToTokens(book):
    if book is not None:
        eachWord = book.lower().split()
        return eachWord
    else:
        return None
        

def unigramModel(tokens):
    hash_map = {}

    if tokens is not None:
        for element in tokens:
           
            word = element
            word = word.replace(",","").replace(".","").replace("/","").replace("?","").replace("(","").replace(")","").replace("[","").replace("]","").replace("-","").replace("_","").replace("'","").replace("!","").replace(":","").replace(";","")
            
            word=unicodedata.normalize('NFKD', word)
           
            list_word=list(word)
            i=0
            lenn=len(word)
            while i<lenn:
                charr=list_word[i]
                if charr.isalpha():
                   
                    if charr in hash_map:
                        hash_map[charr] = hash_map[charr] + 1
                    else:
                        hash_map[charr] = 1
                i+=1        
        return hash_map
    else:
        return None

def bigramModel(tokens):
    hash_map = {'aa':0,'ab':0,'ac':0,'ad':0,'ae':0,'af':0,'ag':0,'ah':0,'ai':0, 'aj':0, 'ak':0, 'al':0, 'am':0, 'an':0, 'ao':0, 'ap':0, 'aq':0, 'ar':0, 'as':0, 'at':0, 'au':0, 'av':0, 'aw':0, 'ax':0, 'ay':0, 'az':0,
                'ba':0,'bb':0,'bc':0,'bd':0,'be':0,'bf':0,'bg':0,'bh':0,'bi':0, 'bj':0, 'bk':0, 'bl':0, 'bm':0, 'bn':0, 'bo':0, 'bp':0, 'bq':0, 'br':0, 'bs':0, 'bt':0, 'bu':0, 'bv':0, 'bw':0, 'bx':0, 'by':0, 'bz':0,
                'ca':0,'cb':0,'cc':0,'cd':0,'ce':0,'cf':0,'cg':0,'ch':0,'ci':0, 'cj':0, 'ck':0, 'cl':0, 'cm':0, 'cn':0, 'co':0, 'cp':0, 'cq':0, 'cr':0, 'cs':0, 'ct':0, 'cu':0, 'cv':0, 'cw':0, 'cx':0, 'cy':0, 'cz':0,
                'da':0,'db':0,'dc':0,'dd':0,'de':0,'df':0,'dg':0,'dh':0,'di':0, 'dj':0, 'dk':0, 'dl':0, 'dm':0, 'dn':0, 'do':0, 'dp':0, 'dq':0, 'dr':0, 'ds':0, 'dt':0, 'du':0, 'dv':0, 'dw':0, 'dx':0, 'dy':0, 'dz':0,
                'ea':0,'eb':0,'ec':0,'ed':0,'ee':0,'ef':0,'eg':0,'eh':0,'ei':0, 'ej':0, 'ek':0, 'el':0, 'em':0, 'en':0, 'eo':0, 'ep':0, 'eq':0, 'er':0, 'es':0, 'et':0, 'eu':0, 'ev':0, 'ew':0, 'ex':0, 'ey':0, 'ez':0,
                'fa':0,'fb':0,'fc':0,'fd':0,'fe':0,'ff':0,'fg':0,'fh':0,'fi':0, 'fj':0, 'fk':0, 'fl':0, 'fm':0, 'fn':0, 'fo':0, 'fp':0, 'fq':0, 'fr':0, 'fs':0, 'ft':0, 'fu':0, 'fv':0, 'fw':0, 'fx':0, 'fy':0, 'fz':0,
                'ga':0,'gb':0,'gc':0,'gd':0,'ge':0,'gf':0,'gg':0,'gh':0,'gi':0, 'gj':0, 'gk':0, 'gl':0, 'gm':0, 'gn':0, 'go':0, 'gp':0, 'gq':0, 'gr':0, 'gs':0, 'gt':0, 'gu':0, 'gv':0, 'gw':0, 'gx':0, 'gy':0, 'gz':0,
                'ha':0,'hb':0,'hc':0,'hd':0,'he':0,'hf':0,'hg':0,'hh':0,'hi':0, 'hj':0, 'hk':0, 'hl':0, 'hm':0, 'hn':0, 'ho':0, 'hp':0, 'hq':0, 'hr':0, 'hs':0, 'ht':0, 'hu':0, 'hv':0, 'hw':0, 'hx':0, 'hy':0, 'hz':0,
                'ia':0,'ib':0,'ic':0,'id':0,'ie':0,'if':0,'ig':0,'ih':0,'ii':0, 'ij':0, 'ik':0, 'il':0, 'im':0, 'in':0, 'io':0, 'ip':0, 'iq':0, 'ir':0, 'is':0, 'it':0, 'iu':0, 'iv':0, 'iw':0, 'ix':0, 'iy':0, 'iz':0,
                'ja':0,'jb':0,'jc':0,'jd':0,'je':0,'jf':0,'jg':0,'jh':0,'ji':0, 'jj':0, 'jk':0, 'jl':0, 'jm':0, 'jn':0, 'jo':0, 'jp':0, 'jq':0, 'jr':0, 'js':0, 'jt':0, 'ju':0, 'jv':0, 'jw':0, 'jx':0, 'jy':0, 'jz':0,
                'ka':0,'kb':0,'kc':0,'kd':0,'ke':0,'kf':0,'kg':0,'kh':0,'ki':0, 'kj':0, 'kk':0, 'kl':0, 'km':0, 'kn':0, 'ko':0, 'kp':0, 'kq':0, 'kr':0, 'ks':0, 'kt':0, 'ku':0, 'kv':0, 'kw':0, 'kx':0, 'ky':0, 'kz':0,
                'la':0,'lb':0,'lc':0,'ld':0,'le':0,'lf':0,'lg':0,'lh':0,'li':0, 'lj':0, 'lk':0, 'll':0, 'lm':0, 'ln':0, 'lo':0, 'lp':0, 'lq':0, 'lr':0, 'ls':0, 'lt':0, 'lu':0, 'lv':0, 'lw':0, 'lx':0, 'ly':0, 'lz':0,
                'ma':0,'mb':0,'mc':0,'md':0,'me':0,'mf':0,'mg':0,'mh':0,'mi':0, 'mj':0, 'mk':0, 'ml':0, 'mm':0, 'mn':0, 'mo':0, 'mp':0, 'mq':0, 'mr':0, 'ms':0, 'mt':0, 'mu':0, 'mv':0, 'mw':0, 'mx':0, 'my':0, 'mz':0,
                'na':0,'nb':0,'nc':0,'nd':0,'ne':0,'nf':0,'ng':0,'nh':0,'ni':0, 'nj':0, 'nk':0, 'nl':0, 'nm':0, 'nn':0, 'no':0, 'np':0, 'nq':0, 'nr':0, 'ns':0, 'nt':0, 'nu':0, 'nv':0, 'nw':0, 'nx':0, 'ny':0, 'nz':0,
                'oa':0,'ob':0,'oc':0,'od':0,'oe':0,'of':0,'og':0,'oh':0,'oi':0, 'oj':0, 'ok':0, 'ol':0, 'om':0, 'on':0, 'oo':0, 'op':0, 'oq':0, 'or':0, 'os':0, 'ot':0, 'ou':0, 'ov':0, 'ow':0, 'ox':0, 'oy':0, 'oz':0,
                'pa':0,'pb':0,'pc':0,'pd':0,'pe':0,'pf':0,'pg':0,'ph':0,'pi':0, 'pj':0, 'pk':0, 'pl':0, 'pm':0, 'pn':0, 'po':0, 'pp':0, 'pq':0, 'pr':0, 'ps':0, 'pt':0, 'pu':0, 'pv':0, 'pw':0, 'px':0, 'py':0, 'pz':0,
                'qa':0,'qb':0,'qc':0,'qd':0,'qe':0,'qf':0,'qg':0,'qh':0,'qi':0, 'qj':0, 'qk':0, 'ql':0, 'qm':0, 'qn':0, 'qo':0, 'qp':0, 'qq':0, 'qr':0, 'qs':0, 'qt':0, 'qu':0, 'qv':0, 'qw':0, 'qx':0, 'qy':0, 'qz':0,
                'ra':0,'rb':0,'rc':0,'rd':0,'re':0,'rf':0,'rg':0,'rh':0,'ri':0, 'rj':0, 'rk':0, 'rl':0, 'rm':0, 'rn':0, 'ro':0, 'rp':0, 'rq':0, 'rr':0, 'rs':0, 'rt':0, 'ru':0, 'rv':0, 'rw':0, 'rx':0, 'ry':0, 'rz':0,
                'sa':0,'sb':0,'sc':0,'sd':0,'se':0,'sf':0,'sg':0,'sh':0,'si':0, 'sj':0, 'sk':0, 'sl':0, 'sm':0, 'sn':0, 'so':0, 'sp':0, 'sq':0, 'sr':0, 'ss':0, 'st':0, 'su':0, 'sv':0, 'sw':0, 'sx':0, 'sy':0, 'sz':0,
                'ta':0,'tb':0,'tc':0,'td':0,'te':0,'tf':0,'tg':0,'th':0,'ti':0, 'tj':0, 'tk':0, 'tl':0, 'tm':0, 'tn':0, 'to':0, 'tp':0, 'tq':0, 'tr':0, 'ts':0, 'tt':0, 'tu':0, 'tv':0, 'tw':0, 'tx':0, 'ty':0, 'tz':0,
                'ua':0,'ub':0,'uc':0,'ud':0,'ue':0,'uf':0,'ug':0,'uh':0,'ui':0, 'uj':0, 'uk':0, 'ul':0, 'um':0, 'un':0, 'uo':0, 'up':0, 'uq':0, 'ur':0, 'us':0, 'ut':0, 'uu':0, 'uv':0, 'uw':0, 'ux':0, 'uy':0, 'uz':0,
                'va':0,'vb':0,'vc':0,'vd':0,'ve':0,'vf':0,'vg':0,'vh':0,'vi':0, 'vj':0, 'vk':0, 'vl':0, 'vm':0, 'vn':0, 'vo':0, 'vp':0, 'vq':0, 'vr':0, 'vs':0, 'vt':0, 'vu':0, 'vv':0, 'vw':0, 'vx':0, 'vy':0, 'vz':0,
                'wa':0,'wb':0,'wc':0,'wd':0,'we':0,'wf':0,'wg':0,'wh':0,'wi':0, 'wj':0, 'wk':0, 'wl':0, 'wm':0, 'wn':0, 'wo':0, 'wp':0, 'wq':0, 'wr':0, 'ws':0, 'wt':0, 'wu':0, 'wv':0, 'ww':0, 'wx':0, 'wy':0, 'wz':0,
                'xa':0,'xb':0,'xc':0,'xd':0,'xe':0,'xf':0,'xg':0,'xh':0,'xi':0, 'xj':0, 'xk':0, 'xl':0, 'xm':0, 'xn':0, 'xo':0, 'xp':0, 'xq':0, 'xr':0, 'xs':0, 'xt':0, 'xu':0, 'xv':0, 'xw':0, 'xx':0, 'xy':0, 'xz':0,
                'ya':0,'yb':0,'yc':0,'yd':0,'ye':0,'yf':0,'yg':0,'yh':0,'yi':0, 'yj':0, 'yk':0, 'yl':0, 'ym':0, 'yn':0, 'yo':0, 'yp':0, 'yq':0, 'yr':0, 'ys':0, 'yt':0, 'yu':0, 'yv':0, 'yw':0, 'yx':0, 'yy':0, 'yz':0,
                'za':0,'zb':0,'zc':0,'zd':0,'ze':0,'zf':0,'zg':0,'zh':0,'zi':0, 'zj':0, 'zk':0, 'zl':0, 'zm':0, 'zn':0, 'zo':0, 'zp':0, 'zq':0, 'zr':0, 'zs':0, 'zt':0, 'zu':0, 'zv':0, 'zw':0, 'zx':0, 'zy':0, 'zz':0}
    
    print(hash_map)

    if tokens is not None:
        for element in tokens:
           
            word = element
            word = word.replace(",","").replace(".","").replace("/","").replace("?","").replace("(","").replace(")","").replace("[","").replace("]","").replace("-","").replace("_","").replace("'","").replace("!","").replace(":","").replace(";","")
            word=unicodedata.normalize('NFKD', word)
            
            list_word=list(word)
            i=0
            lenn=len(word)
            if lenn==1:
                if list_word[0].isalpha():
                    charr="<s>"+list_word[0]
                    if charr in hash_map:
                        hash_map[charr] = hash_map[charr] + 1
                        singleCount['<s>'] +=1
                    else:
                        hash_map[charr] = 1
                        singleCount['<s>'] = 1
            else:                
                while i<lenn-1:
                    charr=list_word[i]+list_word[i+1]
                    if charr.isalpha():
                   
                        if charr in hash_map:
                            hash_map[charr] = hash_map[charr] + 1
                        else:
                            hash_map[charr] = 1
                            
                        if list_word[i] in singleCount:
                            singleCount[list_word[i]] = singleCount[list_word[i]] + 1
                        else:
                            singleCount[list_word[i]] = 1
                    i+=1        
        return hash_map
    else:
        return None


#Unigram Model - English Dump

file = open('en-moby-dick.txt', 'r')
book = file.read()
returnedWords = convertToTokens(book)
ENmap1 = unigramModel(returnedWords)

#print(ENmap1)
count1=sum(ENmap1[item] for item in ENmap1)
#print(count1)


file = open('en-the-little-prince.txt', 'r')
book = file.read()
words = book.lower().split()
ENmap2=unigramModel(words)

#print(ENmap2)
count2=sum(ENmap2[item] for item in ENmap2)
#print(count2)

ENmap={}

A=Counter(ENmap1)
B=Counter(ENmap2)
ENmap=A+B

#print(ENmap)

ENprob=dict()

for k,v in ENmap.items():
    ENprob[k]=(ENmap.get(k)+0.5)/((count1+count2)+(0.5*26))
    
#print(ENprob)


fout = "unigramEN.txt"
fo = open(fout, "w")

for k, v in sorted(ENprob.items()):
    fo.write("("+str(k)+")" + ' = '+ str(v) + '\n')

fo.close()

print("\n")
print("\n")



#Unigram Model - French Dump

file = open('fr-le-petit-prince.txt', 'r')
book = file.read()
words = convertToTokens(book)
FRmap1 = unigramModel(words)

#print(FRmap1)
count1=sum(FRmap1[item] for item in FRmap1)
#print(count1)


file = open('fr-vingt-mille-lieues-sous-les-mers.txt', 'r')
book = file.read()
words = book.lower().split()
FRmap2=unigramModel(words)

#print(FRmap2)
count2=sum(FRmap2[item] for item in FRmap2)
#print(count2)

FRmap={}

A=Counter(FRmap1)
B=Counter(FRmap2)
FRmap=A+B

#print(FRmap)

FRprob=dict()

for k,v in FRmap.items():
    FRprob[k]=(FRmap.get(k)+0.5)/((count1+count2)+(0.5*26))
    
#print(FRprob)


fout = "unigramFR.txt"
fo = open(fout, "w")

for k, v in sorted(FRprob.items()):
    fo.write("("+str(k)+")" + ' = '+ str(v) + '\n')

fo.close()

#Unigram Model - Third Language Dump

file = open('dutch-ebook.txt', 'r', encoding="utf8")
book = file.read()
words = convertToTokens(book)
OTmap1 = unigramModel(words)

#print(OTmap1)
count1=sum(OTmap1[item] for item in OTmap1)
#print(count1)

OTprob=dict()

for k,v in OTmap1.items():
    OTprob[k]=(OTmap1.get(k)+0.5)/((count1)+(0.5*26))
    
#print(OTprob)

fout = "unigramOT.txt"
fo = open(fout, "w")

for k, v in sorted(OTprob.items()):
    fo.write("("+str(k)+")" + ' = '+ str(v) + '\n')

fo.close()







#Bigram Model - English

file = open('en-moby-dick.txt', 'r')
book = file.read()
words = convertToTokens(book)
ENmap1_bi = bigramModel(words)

print(ENmap1_bi)
count1=sum(ENmap1_bi[item] for item in ENmap1_bi)
#print(count1)

file = open('en-the-little-prince.txt', 'r')
book = file.read()
words = book.lower().split()
ENmap2_bi=bigramModel(words)

print(ENmap2_bi)
count2=sum(ENmap2_bi[item] for item in ENmap2_bi)
#print(count2)

ENmap_bi={}

A=Counter(ENmap1_bi)
B=Counter(ENmap2_bi)
ENmap_bi=A+B

print(ENmap_bi)


ENprob_bi=dict()

for k,v in ENmap_bi.items():
    toCalculate=list(k)
    if len(toCalculate) == 2:
        ENprob_bi[k]=(ENmap_bi.get(k)+0.5)/((singleCount.get(toCalculate[0]))+(len(ENmap_bi.keys())*len(ENmap_bi.keys())*0.5))
    else:
        ENprob_bi[k]=(ENmap_bi.get(k)+0.5)/((singleCount.get("<s>"))+(len(ENmap_bi.keys())*len(ENmap_bi.keys())*0.5))
   
#print(ENprob_bi)

#print(singleCount)

fout = "bigramEN.txt"
fo = open(fout, "w")

for k, v in sorted(ENprob_bi.items()):
    toWrite=list(k)
    if len(toWrite)==2:
        fo.write("("+toWrite[1]+"|"+toWrite[0]+")" + ' = '+ str(v) + '\n')
    else:
        fo.write("("+toWrite[3]+"|"+"<"+toWrite[1]+">"+")" + ' = '+ str(v) + '\n')

fo.close()


#Bigram Model - French

#print('\n')
#print('\n')
#print("Bigram - French")

file = open('fr-le-petit-prince.txt', 'r')
book = file.read()
words = convertToTokens(book)
FRmap1_bi = bigramModel(words)

#print(FRmap1_bi)
count1=sum(FRmap1_bi[item] for item in FRmap1_bi)
#print(count1)

file = open('fr-vingt-mille-lieues-sous-les-mers.txt', 'r')
book = file.read()
words = book.lower().split()
FRmap2_bi=bigramModel(words)

#print(FRmap2_bi)
count2=sum(FRmap2_bi[item] for item in FRmap2_bi)
#print(count2)

FRmap_bi={}

A=Counter(FRmap1_bi)
B=Counter(FRmap2_bi)
FRmap_bi=A+B

print(ENmap_bi)


FRprob_bi=dict()

for k,v in FRmap_bi.items():
    toCalculate=list(k)
    if len(toCalculate) == 2:
        FRprob_bi[k]=(FRmap_bi.get(k)+0.5)/((singleCount.get(toCalculate[0]))+(len(FRmap_bi.keys())*len(FRmap_bi.keys())*0.5))
    else:
        FRprob_bi[k]=(FRmap_bi.get(k)+0.5)/((singleCount.get("<s>"))+(len(FRmap_bi.keys())*len(FRmap_bi.keys())*0.5))
   
#print(FRprob_bi)

#print(singleCount)

fout = "bigramFR.txt"
fo = open(fout, "w")

for k, v in sorted(FRprob_bi.items()):
    toWrite=list(k)
    if len(toWrite)==2:
        fo.write("("+toWrite[1]+"|"+toWrite[0]+")" + ' = '+ str(v) + '\n')
    else:
        fo.write("("+toWrite[3]+"|"+"<"+toWrite[1]+">"+")" + ' = '+ str(v) + '\n')

fo.close()


#Bigram Model - Other

file = open('dutch-ebook.txt', 'r', encoding="utf8")
book = file.read()
words = convertToTokens(book)
OTmap1_bi = bigramModel(words)

#print(OTmap1_bi)
count1=sum(OTmap1_bi[item] for item in OTmap1_bi)
#print(count1)

OTprob_bi=dict()

for k,v in OTmap1_bi.items():
    toCalculate=list(k)
    if len(toCalculate) == 2:
        OTprob_bi[k]=(OTmap1_bi.get(k)+0.5)/((singleCount.get(toCalculate[0]))+(len(OTmap1_bi.keys())*len(OTmap1_bi.keys())*0.5))
    else:
        OTprob_bi[k]=(OTmap1_bi.get(k)+0.5)/((singleCount.get("<s>"))+(len(OTmap1_bi.keys())*len(OTmap1_bi.keys())*0.5))

fout = "bigramOT.txt"
fo = open(fout, "w")

for k, v in sorted(OTprob_bi.items()):
    toWrite=list(k)
    if len(toWrite)==2:
        fo.write("("+toWrite[1]+"|"+toWrite[0]+")" + ' = '+ str(v) + '\n')
    else:
        fo.write("("+toWrite[3]+"|"+"<"+toWrite[1]+">"+")" + ' = '+ str(v) + '\n')

fo.close()






with open('Input.txt') as f:
    lines = f.readlines()
    i=0
    while(i<len(lines)):
        j=0
        uni_en_prob=0
        uni_fr_prob=0
        uni_ot_prob=0
        bi_en_prob=0
        bi_fr_prob=0
        bi_ot_prob=0
        #uni_en_prob=1
        #uni_fr_prob=1
        fo = open("out"+str(i+1)+".txt", "w")
        fo.write(lines[i])
        print(lines[i])
        fo.write("\n")
        fo.write("UNIGRAM MODEL:")
        fo.write("\n")
        fo.write("\n")
        eachWord=lines[i]
        #print(eachWord[0])
        while(j<len(lines[i])):
            eachCharac=eachWord[j]
            if eachCharac.isalpha():
                
                if ENprob.get(eachCharac.lower()) is None:
                    uni_en_prob+=0
                    prob_en_uni=0;
                else:
                    uni_en_prob+=math.log10(ENprob.get(eachCharac.lower()))
                    prob_en_uni=math.log10(ENprob.get(eachCharac.lower()))
                
                if FRprob.get(eachCharac.lower()) is None:
                    uni_fr_prob+=0
                    prob_fr_uni=0;
                else:
                    uni_fr_prob+=math.log10(FRprob.get(eachCharac.lower()))
                    prob_fr_uni=math.log10(FRprob.get(eachCharac.lower()))
                
                if OTprob.get(eachCharac.lower()) is None:
                    uni_ot_prob+=0
                    prob_ot_uni=0;
                else:
                    uni_ot_prob+=math.log10(OTprob.get(eachCharac.lower()))
                    prob_ot_uni=math.log10(OTprob.get(eachCharac.lower()))
                    
                fo.write("UNIGRAM: "+eachCharac.lower())
                fo.write("\n")
                fo.write("FRENCH: P("+eachCharac.lower()+") = "+str(prob_fr_uni)+"  ==> log prob of sentence so far: "+str(uni_fr_prob))
                fo.write("\n")
                fo.write("ENGLISH: P("+eachCharac.lower()+") = "+str(prob_en_uni)+"  ==> log prob of sentence so far: "+str(uni_en_prob))
                fo.write("\n")
                fo.write("OTHER: P("+eachCharac.lower()+") = "+str(prob_ot_uni)+"  ==> log prob of sentence so far: "+str(uni_ot_prob))
                fo.write("\n")
                fo.write("\n")
                j+=1
            else:
                j+=1
        j=0
        if uni_en_prob>uni_fr_prob and uni_en_prob>uni_ot_prob:
            fo.write("According to the unigram model, the sentence is in English")
            fo.write("\n")
            print("According to the unigram model, the sentence is in English")
            fo.write("----------------")
        elif uni_fr_prob>uni_en_prob and uni_fr_prob>uni_ot_prob:
            fo.write("According to the unigram model, the sentence is in French")
            fo.write("\n")
            print("According to the unigram model, the sentence is in French")
            fo.write("----------------")
        else:
            fo.write("According to the unigram model, the sentence is in Dutch")
            fo.write("\n")
            print("According to the unigram model, the sentence is in Dutch")
            fo.write("----------------")
        toString = str(lines[i])
        
        eachStringWord = toString.split(" ")
        
        k=0
        fo.write("\n")
        fo.write("BIGRAM MODEL:")
        fo.write("\n")
        fo.write("\n")
        #print(len(eachStringWord))
        while(k<len(eachStringWord)):
            if len(eachStringWord[k]) == 1:
                bi_word=eachStringWord[k]
                
                if ENprob_bi.get("<s>"+bi_word.lower()) is None:
                    bi_en_prob+=0
                    prob_en=0
                else:
                    bi_en_prob+=math.log10(ENprob_bi.get("<s>"+bi_word.lower()))
                    prob_en=math.log10(ENprob_bi.get("<s>"+bi_word.lower()))
                    
                if FRprob_bi.get("<s>"+bi_word.lower()) is None:
                    bi_fr_prob+=0
                    prob_fr=0
                else:
                    bi_fr_prob+=math.log10(FRprob_bi.get("<s>"+bi_word.lower()))
                    prob_fr=math.log10(FRprob_bi.get("<s>"+bi_word.lower()))
                    
                if OTprob_bi.get("<s>"+bi_word.lower()) is None:
                    bi_ot_prob+=0
                    prob_ot=0
                else:
                    bi_ot_prob+=math.log10(OTprob_bi.get("<s>"+bi_word.lower()))
                    prob_ot=math.log10(OTprob_bi.get("<s>"+bi_word.lower()))      
                
                
                fo.write("BIGRAM: "+bi_word.lower())
                fo.write("\n")
                fo.write("FRENCH: P("+bi_word.lower()+"|"+"<s>) = "+str(prob_fr)+"  ==> log prob of sentence so far: "+str(bi_fr_prob))
                fo.write("\n")
                fo.write("ENGLISH: P("+bi_word.lower()+"|"+"<s>) = "+str(prob_en)+"  ==> log prob of sentence so far: "+str(bi_en_prob))
                fo.write("\n")
                fo.write("OTHER: P("+bi_word.lower()+"|"+"<s>) = "+str(prob_ot)+"  ==> log prob of sentence so far: "+str(bi_ot_prob))
                fo.write("\n")
                fo.write("\n")
                
            else:
                bi_word=eachStringWord[k]
                #print(bi_word)
                bi_word=bi_word.lower()
                bi_word=bi_word.replace(",","").replace(".","").replace("/","").replace("?","").replace("(","").replace(")","").replace("[","").replace("]","").replace("-","").replace("_","").replace("'","").replace("!","").replace(":","").replace(";","") 
                bi_word=bi_word.strip()
                #print(bi_word+" empty")
                while j<len(bi_word)-1:
                    charr=bi_word[j]+bi_word[j+1]
                    if FRprob_bi.get(charr) is None:
                        bi_fr_prob+=0
                        prob_fr=0
                    else:  
                        bi_fr_prob+=math.log10(FRprob_bi.get(charr))  
                        prob_fr=math.log10(FRprob_bi.get(charr))
                        
                    if ENprob_bi.get(charr) is None:
                        bi_en_prob+=0
                        prob_en=0
                    else:
                        bi_en_prob+=math.log10(ENprob_bi.get(charr))
                        prob_en=math.log10(ENprob_bi.get(charr))  
                        
                    if OTprob_bi.get(charr) is None:
                        bi_ot_prob+=0
                        prob_ot=0
                    else:
                        bi_ot_prob+=math.log10(OTprob_bi.get(charr))
                        prob_ot=math.log10(OTprob_bi.get(charr)) 
                                     
                    
                    
                    fo.write("BIGRAM: "+charr)
                    fo.write("\n")
                    fo.write("FRENCH: P("+bi_word[j+1]+"|"+bi_word[j]+") = "+str(prob_fr)+"  ==> log prob of sentence so far: "+str(bi_fr_prob))
                    fo.write("\n")
                    fo.write("ENGLISH: P("+bi_word[j+1]+"|"+bi_word[j]+") = "+str(prob_en)+"  ==> log prob of sentence so far: "+str(bi_en_prob))
                    fo.write("\n")
                    fo.write("OTHER: P("+bi_word[j+1]+"|"+bi_word[j]+") = "+str(prob_ot)+"  ==> log prob of sentence so far: "+str(bi_ot_prob))
                    fo.write("\n")
                    fo.write("\n")
                    j+=1
            k+=1
            j=0    
        
        i+=1
        if bi_en_prob>bi_fr_prob and bi_en_prob>bi_ot_prob:
            fo.write("According to the bigram model, the sentence is in English")
            fo.write("\n")
            print("According to the bigram model, the sentence is in English")
            print('\n')
            
        elif bi_fr_prob>bi_en_prob and bi_fr_prob>bi_ot_prob:
            fo.write("According to the bigram model, the sentence is in French")
            fo.write("\n")
            print("According to the bigram model, the sentence is in French")
            print('\n')
        else:
            fo.write("According to the bigram model, the sentence is in Dutch")
            fo.write("\n")
            print("According to the bigram model, the sentence is in Dutch")
            print('\n')
        
    fo.close()
        
            
            
  
    

    