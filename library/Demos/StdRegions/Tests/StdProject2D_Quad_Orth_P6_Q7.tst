<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject2D Quadrilateral Orthonormal basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s quadrilateral -b Ortho_A Ortho_A -o 6 6 -p 7 7</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">2.85686e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.06581e-14</value>
        </metric>
    </metrics>
</test>


