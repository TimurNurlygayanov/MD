from flask import make_response
from flask import render_template
from flask import redirect, url_for


def valid_login(user, password):
    if user == "Timur" and password == "test":
        return True
    return False


def log_out_user():
    resp = redirect(url_for('login'))
    resp.set_cookie('username', '')
    return resp


def log_the_user_in(user):
    resp = redirect(url_for('index'))
    resp.set_cookie('username', user)
    return resp
