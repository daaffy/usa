# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 18:37:23 2013

@author: gordon
"""

class GETag(object):
    
    def __init__(self,*args):
        self.data = dict('firsttag'=args[0] ,'secondtag'= args[1],         
        'size' = args[2],'data' = args[3])
        
    def GetFirstTag(self):
        return self.data('firsttag')        
    def GetSecondTag(self):
        return self.data('secondtag')
    def GetData(self):
        return self.data('data')
    def GetSize(self):
        return self.data('size')
        