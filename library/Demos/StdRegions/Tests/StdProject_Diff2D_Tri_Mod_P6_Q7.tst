<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject_Diff2D Triangle Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s triangle -b Modified_A Modified_B -o 6 6 -p 7 7 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">2.0808e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">2.74891e-13</value>
        </metric>
    </metrics>
</test>


