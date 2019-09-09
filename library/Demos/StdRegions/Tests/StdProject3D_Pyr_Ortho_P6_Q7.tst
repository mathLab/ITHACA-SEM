<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>StdProject3D Pyramid Orthogonal basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s pyramid -b Ortho_A Ortho_A OrthoPyr_C -o 6 6 6 -p 7 7 6</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">5.32907e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">5.32907e-15</value>
        </metric>
    </metrics>
</test>
