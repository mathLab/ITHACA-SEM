<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project2D Quad Lagrange basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7 -c 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.81771e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">1.06581e-14</value>
        </metric>
    </metrics>
</test>
