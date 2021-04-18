import requests as rq
from bs4 import BeautifulSoup
import re
import csv

BASE_URL = r'https://ned.ipac.caltech.edu/byname?objname={}&hconst=67.8&omegam=0.308&omegav=0.692&wmap=4&corr_z=1'

def main():
    with open('galaxies.csv') as f_in, open('results.csv', 'w', newline='') as f_out:
        reader = csv.reader(f_in)
        next(reader)

        writer = csv.writer(f_out)
        writer.writerow(['name', 'number_of_notes', 'number_of_ref'])

        for line in reader:
            name = line[0]
            try:
                number_of_notes, number_of_ref = parse_object(name)
                writer.writerow([name, number_of_notes, number_of_ref])
            except ValueError:
                writer.writerow([name, 'not', 'found'])


def parse_object(name):
    html = get_html(BASE_URL.format(name))
    soup = get_soup(html)

    number_of_notes = soup.find('a', href="#tabs_7")    # сделать "единообразно" не вышло, сайт чуть криво написан, как мне кажется
    number_of_ref   = soup.find('a', "reference")       # поэтому для других вкладок нужно смотреть заново код страницы :(

    number_of_notes = number_of_notes.text.strip()      # перевод объекта bs4 в текстовую строку
    number_of_ref = number_of_ref.text.strip()

    number_of_notes = re.findall('\d+', number_of_notes)  # поиск чисел в строке
    number_of_ref = re.findall('\d+', number_of_ref)
    return(number_of_notes, number_of_ref)


def get_html(url):
    return rq.get(url).content


def get_soup(html):
    return BeautifulSoup(html, "lxml")


if __name__ == '__main__':
    main()
