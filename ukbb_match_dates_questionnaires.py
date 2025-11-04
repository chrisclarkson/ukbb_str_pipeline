import requests
import bs4
import sys
import pandas as pd

index_url = 'https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id='
input=sys.argv[1]

data=pd.read_csv('~/Downloads/compare_cols3.csv')

def get_data(code):
    code=code.replace('p','')
    url=index_url+str(code)
    print(url)
    response = requests.get(url)
    soup = bs4.BeautifulSoup(response.text)
    # print(soup)
    table = soup.find_all('table', attrs={'summary':'Identification'})[0]
    # print(table)
    tr=table.find_all('tr')[1]
    td=tr.find_all('td')[1]
    # print(td)
    a=td.find_all('a', href=True)
    a=a[len(a)-1]
    return(a)

def open_soupling(href):
    response = requests.get(href)
    soup = bs4.BeautifulSoup(response.text)
    table=soup.find_all('table',attrs={'summary':'List of data-fields'})[0]
    trs=table.find_all('tr')
    out=[]
    for tr in trs:
        if 'Field ID' in tr.text:
            continue
        td=tr.find_all('td')
        out.append([td[0].text,td[1].text])
    out=pd.DataFrame(out)
    return(out)

def retrieve_relevant_fields(data,retrieve_df,when_field):
    print(retrieve_df)
    data_sub=data.loc[data['Category.Field.ID'].isin(retrieve_df[0].astype(int))]
    print(when_field)
    when_field_name=retrieve_df.loc[retrieve_df[0].astype(int)==int(when_field)][1].values[0]
    when_field_name='p'+str(when_field)+'_'+when_field_name.replace(' ','_')
    data_sub['date_col']=when_field_name
    print(data_sub)
    data.loc[data['Category.Field.ID'].isin(retrieve_df[0].astype(int))]=data_sub
    return(data)

with open(input,'r') as codes:
	codes=codes.readlines()
	for code in codes:
		code=code.rstrip()
		code=code.split('\t')[0]
		if sum(data['Category.Field.ID'])==0:
			continue
		href=str(get_data(code))
		href=href.split(' href="')[1].split('">')[0]
		when_field=int(href.replace('label.cgi?id=',''))
		href='https://biobank.ndph.ox.ac.uk/ukb/'+href
		retrieve_df=open_soupling(href)
		data=retrieve_relevant_fields(data,retrieve_df,int(code.rstrip().replace('p','')))

print(data)
data.to_csv('test.csv',index=False)