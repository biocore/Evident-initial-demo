from flask import jsonify

def load_studies():
    studies = {}
    studydata = open('data/studies.txt').readlines()
    studydata = studydata[1:]

    # Writting JS headers
    jslines = "{'demo':{"

    for row in studydata:
        if row.startswith("#") or len(row)==0:
           continue
        vals = row.strip().split("\t")
        jslines += "'%s': ['%s']," % (vals[1],vals[3])
    jslines += "}, 'full':{"

    for row in studydata:
        if row.startswith("#") or len(row)==0:
           continue
        vals = row.strip().split("\t")

        # opening selectors
        selectors = open("data/%s_selectors.txt" % vals[1], "U").read().replace("\n","?")

        jslines += "'%s': ['%s','%s']," % (vals[1],vals[0],selectors)

    jslines = jslines[:-1]
    jslines += '}}'
    return jsonify(result=jslines)

