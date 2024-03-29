# Superconductivity_research　Vortexratcheteffect_simulation


 ## Memo
### 第一マッチングフィールドが一番正味の起電力が高い<br>
→第一マッチングフィールドから第二マッチングフィールド間においてはボルテックスの数密度の増加によって起電力そのものは増加するがRachetBになって正味の起電力は低くなる。<br>
第一マッチングフィールド以下においてはボルテックスの平均速度は変化しないがボルテックスの数密度が減少するので起電力は減少する<br>

### ピーク電流（最大の起電力が期待できる電流値）においてサイト間距離は短ければ短いほど期待できる正味の起電力は大きくなる<br>
→サイト間距離が延びるとボルテックスの相互作用が弱まることで平均速度が遅くなる。かつボルテックスの数密度が減少するため<br>

### サイト間距離が長くなるとオンセット電流（ピニングサイトが外れる最小の電流値）も増加する（臨界電流に近づく）<br>
→サイト間距離が長くなるとボルテックスの相互作用が弱まり、ピニングサイトを外すために必要なローレンツ力が増加するため。<br>

### ラチェット効果を起こす重要な要素<br>
-非対称性の作り方（形の違い）<br>
-配列（連続するピニングサイトの距離や形状）<br>

### 大中小円より大小円２つのピニングサイトで非対称性を作った場合がラチェト効果の観測により有効<br>
→triple最小値が-0.312、Double の最小値が-0.385でDouble のほうが平均速度が強い<br>

## 調べる値
①	ピニングサイト位置<br>
グラフ横軸FL/Fp：ローレンツ力とピニング力の比VS縦軸sitedistance、平均速度を色でxy座標に盛り込む。（カンターマップで３つの変数で表す）<br>


②	ピニングサイトの大きさ、サイト間距離依存<br>
大円RL＝１～３<br>
小円Rs＝１<br>
Sitedistance=0~5<br>
DpY=16<br>
(1＝100nm換算)<br>

グラフ横軸FL/Fp：ローレンツ力とピニング力の比<br>
縦軸sitedistance、平均速度を色でxy座標に盛り込む。（カンターマップで３つの変数で表す）<br>

③	電流値・起電力の換算<br>
②の条件の結果において各条件における平均速度のオンセットの位置におけるオンセット電流と最大平均速度のピーク電流を臨界電流と比の形で求める。<br>
最大平均速度から正味の起電力を求める<br>

グラフ：横軸sitedistance,縦軸：J/Jc電流値と臨界電流の比、ピーク電流とオンセット電流をマッピング<br>
グラフ：横軸sitedistance,縦軸：Enet(正味の起電力)、ピーク電流における正味の起電力<br>

④	マッチングフィールド（磁束密度の値、ボルテックスを一つずつ増減させることでマッチングフィールドからの磁場の割合的なずれを考えるため）<br>
グラフ：横軸FL/Fp：ローレンツ力とピニング力の比、縦軸：H/H1第一マッチングフィールドと磁場の比、平均速度を色でxy座標に盛り込む（カンターマップ）<br>
グラフ：横軸：H/H1第一マッチングフィールドと磁場の比、縦軸：Enet正味の起電力<br>

## 実験より分かった条件
①	より<br>
対称性がラチェット効果に影響を及ぼすと分かった<br>
Rl＝2,Rs=1,siteDistance=3～9,DpY=16でsitedistance=6のときボルテックスの平均速度反転<br>

②	より<br>
Rl＝３Rs＝１が一番顕著に平均速度に違いがある（Rl＝１～３、Rs＝３）Sitedistance=0.2<br>

③	より<br>
Rl＝３,Rs=1,sitedistance=0～5で変化<br>
Rl=3,Rs=1,J/Jc=0.977,sitedistance=0.2で最大平均速度<br>
Rl=3,Rs=1,J/Jc=0.977～0.996,sitedistance=0.2~1.7で最大平均速度は上昇（振動してるけど）<br>
Rl=3,Rs=1,J/Jc=0.996,sitedistance=1.7～5.0で平均速度ほぼ変わらない<br>
Rl＝３,Rs=1,sitedistance=0.2で最大の正味の起電力E(net)=6.78(mV)(sitedistance=0～5)<br>

④	より<br>
H/H1=0.5～1.25付近,Fl/Fp=0.970～0.992付近で負のラチェット効果発現、H/H1＝１で最大正味の起電力0.053(V)(H/H1=0.5～2)<br>
