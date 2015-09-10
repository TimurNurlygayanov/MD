import json
from flask import Flask
from flask import jsonify
from flask import render_template
from flask import request
from flask import redirect, url_for
from flask_bootstrap import Bootstrap

import tools


app = Flask(__name__)
Bootstrap(app)


def check_auth(request):
    username = request.cookies.get('username')
    if username:
        return render_template('index.html', username=username)
    else:
        return redirect(url_for('login'))


@app.route('/experiment_list', methods=['GET', 'POST'])
def experiment_list():
    #check_auth(request)
    result = {'total': 200,
              'rows': [ {'name': 'test', 'author': 'test',
                         'date': 'Apr 1 1988', 'description': 'test'},
                        {'name': 'test2', 'author': 'test',
                         'date': 'Apr 2 1988', 'description': 'test2'}
                      ]
             }
    #return json.dumps(result)
    return jsonify(result)


@app.route("/")
def index():
    return check_auth(request)


@app.route('/logout', methods=['GET'])
def logout():
    return tools.log_out_user()


@app.route('/login', methods=['POST', 'GET'])
def login():
    error = None
    if request.method == 'GET':
        return render_template('login.html')
    if request.method == 'POST':
        if tools.valid_login(request.form.get('username'),
                             request.form.get('password')):
           return tools.log_the_user_in(request.form['username'])
        else:
            error = 'Invalid username/password'
    return render_template('login.html', error=error)


if __name__ == "__main__":
    app.run(host='0.0.0.0', debug=True)
