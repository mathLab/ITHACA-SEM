<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject_Diff2D Quadrilateral Lagrange basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">4.14376e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.13687e-13</value>
        </metric>
    </metrics>
</test>


