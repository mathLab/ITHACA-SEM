<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject2D Quadrilateral Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s quadrilateral -b Modified_A Modified_A -o 6 6 -p 7 7</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.08794e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">8.66862e-13</value>
        </metric>
    </metrics>
</test>


