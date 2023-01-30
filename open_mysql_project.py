#!/usr/bin/env python
# -*- coding:utf-8 -*-

import mysql.connector
from urllib.parse import urlparse


def openproject():
    urlsql = urlparse('mysql://root:kwmjbqb9py@localhost:3306/DIBproject')

    conn = mysql.connector.connect(
        host=urlsql.hostname or 'localhost',
        port=urlsql.port or 3306,
        user=urlsql.username or 'root',
        password=urlsql.password or 'kwmjbqb9py',
        database=urlsql.path[1:],
    )

    cur = conn.cursor()

    return conn, cur
